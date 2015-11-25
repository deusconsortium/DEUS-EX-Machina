#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  sans titre.py
#
#  Copyright 2015 pdefromont <pdefromont@advil.obspm.fr>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#


from DEUSGraphics import *

test = DEUSGraphics(1100.)
# if working in PyCharm
test.setDataPath("/home/pdefromont/advil/data/home/pdefromont/code/DEUS-EX-Machina/ProfileTracer/data/output/")
root_file = "boxlen2592_npart2048_lcdmw5_output"

outputs = [56, 52, 45, 39, 25, 18, 6, 1]
rz = num.zeros(size(outputs))
dz = num.zeros(size(outputs))
a = num.zeros(size(outputs))

r1min = 5.0
r1max = 25.0

for i in range(size(outputs)):
    file_name = root_file + str(outputs[i]) + "_Zwork"
    test.Load(file_name)
    a[i] = test._a

    figure(2)
    grid(True)
    r, mf, sf, mv, sv = test.getMeanProfile(R1=15.)
    plot(r, mf, color=rbcolor(i, size(outputs)), label="z = " + str(1. / a[i] - 1.))

    if i == 0:
        N = test.getProfileNumber()
        r10 = test._r1_full
        tab = num.ones(N)
        rz[0] = 1.0
        dz[0] = 0.0
    else:
        r1 = test._r1_full
        for j in range(size(r1)):
            if r1[j] is not None and r10[j] is not None and r1max > r1[j] > r1min and r1max > r10[j] > r1min:
                tab[j] = r1[j] / r10[j]
            else:
                tab[j] = None
        figure(0)
        grid(True)
        bins = num.linspace(0.0, 5.0, 41)
        hist(tab[tab != num.array(None)], bins, alpha=0.3, color=rbcolor(i, size(outputs)), normed=True,
             histtype="stepfilled", label="z = " + str(1.0 / a[i] - 1.0))
        legend()
        xlabel("$R_1(z)/R_1(0)$")
        ylabel("$P$")
        rz[i] = num.nanmean(tab)
        dz[i] = num.nanstd(tab)

figure(1)
grid(True)
ylabel("$R_1(a)/R_1(0)$")
xlabel("$a$")
errorplot(a, rz, dz)
ylim((0.0, 1.5))
figure(2)
ylabel("$\\Delta(r) + 1$")
xlabel("$r$ in $[Mpc.h^{-1}]$")
legend()
show()
