#define _XOPEN_SOURCE 500
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>

//#define OFFSET_MAX 274877906944 // 4*4096^3
//#define OFFSET_MAX 2199023255552 // 4*8192^3
#define OFFSET_MAX 140737488355328 // 4*32768^3

int parallel_write(int,long,long long,void*);
int parallel_read(int,long,long long,void*);
#if defined ADD0US
void f77_parallel_openread(char *,int *,int *);
void f77_parallel_openwrite(char *,int *,int *);
void f77_parallel_close(int *);
void f77_parallel_read(int*,long*,long long*,void*);
void f77_parallel_write(int*,long*,long long*,void*);
#elif defined ADD2US
void f77_parallel_openread__(char *,int *,int *);
void f77_parallel_openwrite__(char *,int *,int *);
void f77_parallel_close__(int *);
void f77_parallel_read__(int*,long*,long long*,void*);
void f77_parallel_write__(int*,long*,long long*,void*);
#else
void f77_parallel_openread_(char *,int *,int *);
void f77_parallel_openwrite_(char *,int *,int *);
void f77_parallel_close_(int *);
void f77_parallel_read_(int*,long*,long long*,void*);
void f77_parallel_write_(int*,long*,long long*,void*);
#endif


#if defined ADD0US
void f77_parallel_openread(char *filename, int *fnamelen, int *fd) {
#elif defined ADD2US
void f77_parallel_openread__(char *filename, int *fnamelen, int *fd) {
#else
void f77_parallel_openread_(char *filename, int *fnamelen, int *fd) {
#endif
  filename[*fnamelen]='\0';
  if( (*fd = open64(filename,O_RDONLY,0644)) == -1) perror("open");
}

#if defined ADD0US
void f77_parallel_openwrite(char *filename, int *fnamelen, int *fd) {
#elif defined ADD2US
void f77_parallel_openwrite__(char *filename, int *fnamelen, int *fd) {
#else
void f77_parallel_openwrite_(char *filename, int *fnamelen, int *fd) {
#endif
  filename[*fnamelen]='\0';
  if( (*fd = open64(filename,O_RDWR|O_CREAT,S_IRUSR|S_IWUSR)) == -1) perror("open");
}

#if defined ADD0US
void f77_parallel_close(int *fd) {
#elif defined ADD2US
void f77_parallel_close__(int *fd) {
#else
void f77_parallel_close_(int *fd) {
#endif
  if ( close(*fd) == -1 ) perror("close");
}

// This routine writes size bytes of buffer at position offset
// in the file filename.
// Can be used in multiple concurrent access

int parallel_write(int fd, long size, long long offset, void *buffer) {

  ssize_t stat; // I/O status
  off64_t stat2;

  if (size + offset > OFFSET_MAX) {
    fprintf(stderr,"You are trying to access a file location\n");
    fprintf(stderr,"which is bigger than %ld\n",(long long)OFFSET_MAX);
    fprintf(stderr,"Verify your code and/or change OFFSET_MAX\n");
    return(1);
  }

  if( (stat2 = lseek64(fd,offset,SEEK_SET))                  != offset ) perror("lseek");
  if( (stat = write(fd,buffer,size))                         != size )   perror("read");

  if (stat != size) return(1);
  else return(0);
}

int parallel_read(int fd, long size, long long offset, void *buffer) {

  ssize_t stat; // I/O status
  off64_t stat2;

  if (size + offset > OFFSET_MAX) {
    fprintf(stderr,"You are trying to access a file location\n");
    fprintf(stderr,"which is bigger than %ld\n",(long long)OFFSET_MAX);
    fprintf(stderr,"Verify your code and/or change OFFSET_MAX\n");
    return(1);
  }

  if( (stat2 = lseek64(fd,offset,SEEK_SET)) != offset ) perror("lseek");
  if( (stat = read(fd,buffer,size))         != size )   perror("read");

  //printf("read: fd=%d size=%ld offset=%lld stat=%ld stat2=%lld\n",fd,size,offset,stat,stat2);

  if (stat != size) return(1);
  else return(0);
}

// Fortran wrappers

#if defined ADD0US
void f77_parallel_read(int *fd, long *size, long long *offset, void *buffer) {
#elif defined ADD2US
void f77_parallel_read__(int *fd, long *size, long long *offset, void *buffer) {
#else
void f77_parallel_read_(int *fd, long *size, long long *offset, void *buffer) {
#endif
  int stat;
  stat = parallel_read(*fd,*size,*offset,buffer);
  if (stat !=0) fprintf(stderr,"Failure to read with parallel_read\n");
}

#if defined ADD0US
void f77_parallel_write(int *fd, long *size, long long *offset, void *buffer) {
#elif defined ADD2US
void f77_parallel_write__(int *fd, long *size, long long *offset, void *buffer) {
#else
void f77_parallel_write_(int *fd, long *size, long long *offset, void *buffer) {
#endif
  int stat;
  stat = parallel_write(*fd,*size,*offset,buffer);
  if (stat !=0) fprintf(stderr,"Failure to write with parallel_read\n");
}



///
/*
TMP_main(int argc, char **argv) {

  int fd;
  char filename[128];
  long long noffset;
  long nsize;
  float *buffer;
  int i, stat;
  int myid, nproc;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if (argc != 3 && myid==0) {
    fprintf(stderr,"io filename nsize\n");
    MPI_Finalize();
    exit(1);
  }
    

  strcpy(filename,argv[1]);
  nsize=atoi(argv[2]); // number of floats to read/write
  // noffset=atoi(argv[3]); // offset in sizeof(float)
  noffset=myid*nsize;
  
  
  fd = open(filename,O_RDWR|O_CREAT,S_IRUSR|S_IWUSR);
  buffer = (float*)malloc(nsize*sizeof(float));
  for (i=0;i<nsize;i++) buffer[i] = (float)(i+nsize*myid);
  stat=pwrite64(fd,(void*)buffer,nsize*sizeof(float),noffset*sizeof(float));
  close(fd);
  fprintf(stderr,"coucou from proc # %d: stat = %d\n",myid,stat);
  MPI_Barrier(MPI_COMM_WORLD);
  return;
  

  
  fd = open(filename,O_RDONLY);
  buffer = (float*)malloc(nsize*sizeof(float));
  stat=pread64(fd,(void*)buffer,nsize*sizeof(float),noffset*sizeof(float));
  for (i=0;i<nsize;i++) fprintf(stderr,"%d %f\n",i+nsize*myid,buffer[i]);
  close(fd);
  return;
  
  MPI_Finalize();

}
*/
