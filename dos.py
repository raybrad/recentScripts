#***********************************************************
# * Copyright (C) 2013 Alexey V. Akimov
# * This file is distributed under the terms of the
# * GNU General Public License as published by the
# * Free Software Foundation; either version 3 of the
# * License, or (at your option) any later version.
# * http://www.gnu.org/copyleft/gpl.txt
# 
# * Please modify the sections below if you add new functionality
# * to this module (add new "Last modified" entry for each unique
# * contributor, move the "Last modified" entry to the top for
# * each latest modification)
#
# * Initial version: 12/9/2013      by Alexey V. Akimov
# * Last modified: 12/9/2013        by Alexey V. Akimov
# 
#***********************************************************/

import math
import sys

def gaussian_dos(filename,outfile,specs,dos_type,emin,emax,dx,var,spin_method):
# filename - name of the input file (output of Gaussian)
#  Example: "CD28S38_PM6.LOG"
#
# outfile - name of the file into which (p)DOS will be written
#  Example: "Cd28S38_pm6.txt"
#
# specs - species on which to project states
#  Example:  [range(0,28),range(28,66)] - first set of atomic species includes atoms from 0 to 27 (not including 28)
#                                         second set of atomic specied includes atoms from 28 to 65 (not including 66)
#
# dos_type - how to interpred data in specs
#  Possible options:  "atomic" - the indices in the specs at atomic indices
#
# emin, emax - range of the energies for which (p)DOS will be computed (in eV)
#
# dx - (p)DOS grid spacing (in eV)
#
# var - Gaussian broadening width of the DOS lines (in eV)
#
# spin_method: 
#  Possible options:
#    "restricted" - use it if restricted (closed-shell) calculations have been performed
#    "unrestricted" - use it if unrestricted (open-shell) calculations have been performed


    key1 = "Alphaocc.eigenvalues--"
    key2 = "Alphavirt.eigenvalues--"
    key3 = "Betaocc.eigenvalues--"
    key4 = "Betavirt.eigenvalues--"
    key5 = "Grossorbitalpopulations:"  # 3 tokens
    key6 = "TotalAlphaBetaSpin"  # 4 tokens
    key7 = "Condensedtoatoms(allelectrons):"  # 5 tokens
    key8 = "NBasis="  # 1 token
    key9 = "Eigenvalues--"  # 2 tokens

    E_occ = []
    E_alp = []
    At_indx = []
    At_label = []
    Orb_label = []
    Orb_pop = []

    f = open(filename,"r")
    A = f.readlines()   
    f.close()



    #================ Read in all necessary data from Gaussian output file ================
    nit = 0
    if spin_method=="restricted":
        nit = 0
    elif spin_method=="unrestricted":
        nit = 3

    # First determine the number of basis functions (# of states)
    Nbas = 0
    i = 0
    status1 = 0
    while status1==0:
        tmp = A[i].split()
        sz = len(tmp)
        if sz>=2:
            if(tmp[0]==key8):
                Nbas = int(float(tmp[1]))
                print "Number of states = ",Nbas
                status1 = 1
        i = i + 1

    #Prepare storage for MO expansions
    MO = []
    for i in range(0,Nbas):
        mo = []
        for j in range(0,Nbas):
            mo.append(0.0)
        MO.append(mo)



    # Now read the MO expansions
    status1 = 0    
    i = 0  # MO index

    nlines = len(A)
    k = 0
    while k<nlines:
        tmp = A[k].split()
        sz = len(tmp)

        # Reading MO expansion
        if sz>2 and i<Nbas:
            key = tmp[0]+tmp[1]
            if key==key9:
                norbs = sz - 2  # how many orbitals are printed on this line

                k = k+1  # advance to next line

                # Read next Nbas lines
                for k1 in range(0,Nbas):                   
                    tmp = A[k+k1].split()
                    sz = len(tmp)
                    npref = sz - norbs  # how many slots are for prefix

#                    print i,tmp[npref:]

                    for j in range(0,norbs):
#                        print i,j,tmp[npref+j]
                        MO[i+j][k1] = float(tmp[npref+j])
   
            
                i = i + norbs
                k = k + Nbas  # advance to Nbas lines

        

        # Reading energies
        if sz>4:
            key = tmp[0]+tmp[1]+tmp[2]+tmp[3]
            if key==key1:
                for e in tmp[4:]:
                    E_occ.append(float(e)/0.036749309)
            if key==key1 or key==key2:
                for e in tmp[4:]:
                    E_alp.append(float(e)/0.036749309)  # convert from Ha to eV


        # Reading atomic and orbital labels
        if sz==3:
            key = tmp[0]+tmp[1]+tmp[2]
            if key==key5:
                k = k + 2  # skip one line

                at_ind = 1
                at_typ = ""
                pop = 0.0
                orb_lab = ""


                for k1 in range (0,Nbas):
                    tmp = A[k+k1].split()
                    sz = len(tmp)

                    if sz==(5+nit):
                        at_ind = int(float(tmp[1]))
                        at_typ = tmp[2]
                        orb_lab= tmp[3]
                        pop = float(tmp[4])

                    elif sz==(3+nit):
                        orb_lab = tmp[1]
                        pop = float(tmp[2])          

                    elif sz==(4+nit):  # D 0 orbital
                        orb_lab = tmp[1]+" "+tmp[2]
                        pop = float(tmp[3])          

                    At_indx.append(at_ind)
                    At_label.append(at_typ)
                    Orb_label.append(orb_lab)
                    Orb_pop.append(pop)

                k = k + Nbas  # advance k by Nbas lines

        k = k + 1


# Optional - comment if notinterested
    print len(E_alp), E_alp
    print len(At_indx),At_indx
    print len(At_label),At_label
    print len(Orb_label),Orb_label
    print len(Orb_pop),Orb_pop


    #============== Prepare grid ======================================
    Norb = len(E_alp)  # original number of energy levels
    Nocc = len(E_occ)

    E_homo = E_occ[Nocc-1]

    ngrid = int((E_alp[Norb-1] - E_alp[0])/dx)
    npdos = len(specs) # now many pDOSs we want to compute

    X = []
    pDOS = []
    i = 0
    while i<ngrid:
        X.append(E_alp[0] + i * dx)
        
        pdos = []
        j = 0
        while j<npdos:
            pdos.append(0.0)
            j = j + 1
        pDOS.append(pdos)

        i = i + 1

    # so now pDOS[i][j] - contains DOS for j-th projection group (group of atoms perhaps) at i-th grid point
    #============== Now we are ready to make up pDOS ==================
    if(dos_type=="atomic"):

        f1 = open(outfile,"w")

        prefac = 1.0/(var*math.sqrt(2.0*math.pi))
        alp = 0.5/(var*var)

        integ = 0.0 # integral of DOS

        i = 0
        while i<ngrid:

            if emin<=((X[i]-E_homo)) and ((X[i]-E_homo)<=emax):

                j = 0
                line = str(X[i]-E_homo)+"  "
                while j<npdos:
        
                    pDOS[i][j] = 0.0
        
                    #==============================================================
                    ii = 0
                    while ii<Norb:
        
                        wi = 0.0
                        jj=0
                        while jj<Norb:
                            if(At_indx[jj] in specs[j]):
                                wi = wi + MO[ii][jj]**2 # weight of jj-th atomic orbital in ii-th MO
                            jj = jj + 1
        
                        w = prefac*math.exp(-alp*(X[i] - E_alp[ii])**2)  # normalized Gaussian                        
        
                        pDOS[i][j] = pDOS[i][j] + w*wi
        
        
                        ii = ii + 1
                    #==============================================================
                    integ = integ + pDOS[i][j]*dx
        
                    line = line + str(pDOS[i][j])+"  "
                    j = j + 1
                line = line + "\n" 
                f1.write(line)

            i = i + 1

        print "Integral of DOS = ",integ

        f1.close()     


def yaehmop_dos(filename,outfile,specs,dos_type,emin,emax,dx,var,spin_method):
# filename - name of the input file (output of Gaussian)
#  Example: "input_cd28s38_v5.bind.out"
#
# outfile - name of the file into which (p)DOS will be written
#  Example: "Cd28S38_eht_v5.txt"
#
# specs - species on which to project states
#  Example:  [range(0,28),range(28,66)] - first set of atomic species includes atoms from 0 to 27 (not including 28)
#                                         second set of atomic specied includes atoms from 28 to 65 (not including 66)
#
# dos_type - how to interpred data in specs
#  Possible options:  "atomic" - the indices in the specs at atomic indices
#
# emin, emax - range of the energies for which (p)DOS will be computed (in eV)
#
# dx - (p)DOS grid spacing (in eV)
#
# var - Gaussian broadening width of the DOS lines (in eV)
#
# spin_method: 
#  Possible options:
#    "restricted" - use it if restricted (closed-shell) calculations have been performed
#    "unrestricted" - use it if unrestricted (open-shell) calculations have been performed


    key1 = "#Num_Orbitals:"
    key2 = ";***Wavefunctions***" # 4 tokens
    key3 = "#*******Energies(ineV)andOccupationNumbers*******" # 9 tokens
    key4 = "[0.000"  # 1 token

    E_occ = []
    E_alp = []
    At_indx = []
    At_label = []
    Orb_label = []
    Orb_pop = []

    f = open(filename,"r")
    A = f.readlines()   
    f.close()
    nlines = len(A)



    #================ Read in all necessary data from Gaussian output file ================
    nit = 0
    if spin_method=="restricted":
        nit = 0
    elif spin_method=="unrestricted":
        nit = 3

    # First determine the number of basis functions (# of states)
    Nbas = 0
    i = 0
    status1 = 0
    while status1==0 and i<nlines:
        tmp = A[i].split()
        sz = len(tmp)
        if sz>=2:
            if(tmp[0]==key1):
                Nbas = int(float(tmp[1]))
                print "Number of states = ",Nbas
                status1 = 1
        i = i + 1

    #Prepare storage for MO expansions
    MO = []
    for i in range(0,Nbas):
        mo = []
        for j in range(0,Nbas):
            mo.append(0.0)
        MO.append(mo)


    # Now read the MO expansions
    status1 = 0    
    i = 0  # MO index
    k = 0
    while k<nlines:
        tmp = A[k].split()
        sz = len(tmp)

        # Reading MO expansion
        if status1==1:
            norbs = len(tmp)/3  # works as long as index of atom is <100 !!!!!!!
            
            for j in range(0,norbs):
                At_label.append(tmp[3*j+0][:-1])
                At_indx.append(int(float(tmp[3*j+1][:-1])))
                Orb_label.append(tmp[3*j+2])

            k = k + 1
            # Read next Nbas lines
            for k1 in range(0,Nbas):                   
                tmp = A[k+k1].split()

                sz = len(tmp) # must be equal to norbs + 1
                if sz!=(norbs+1):
                    print "# of atoms is >=100 - check the code"
                    sys.exit(0)

                for j in range(0,norbs):
                    MO[k1][i+j] = float(tmp[1+j])
            
            i = i + norbs
            k = k + Nbas-1  # advance to Nbas lines

            if i>=Nbas:
                status1 = 0

        elif status1==0:

            if sz>=4: 
                key = tmp[0]+tmp[1]+tmp[2]+tmp[3]
                if key==key2:
                    status1 = 1   # this shoud be after if status1==1 loop
                    k = k + 1 # to skip one line
                       
            # Reading energies
            if sz==9:
                key = tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+tmp[6]+tmp[7]+tmp[8]
                if key==key3:
                    k = k + 1
                    for k1 in range(0,Nbas):
                        tmp = A[k+k1].split()

                        e = float(tmp[1])
                        if tmp[2]!=key4:
                            E_occ.append(e)
                        E_alp.append(e) 

                    k = k + Nbas  # advance k by Nbas lines
        k = k + 1


# Optional - comment if notinterested
    print len(E_alp), E_alp
    print len(At_indx),At_indx
    print len(At_label),At_label
    print len(Orb_label),Orb_label


    #============== Prepare grid ======================================
    Norb = len(E_alp)  # original number of energy levels
    Nocc = len(E_occ)

    E_homo = E_occ[Nocc-1]

    ngrid = int((E_alp[Norb-1] - E_alp[0])/dx)
    npdos = len(specs) # now many pDOSs we want to compute

    print "Norb= ",Norb," Nocc= ",Nocc, " E_homo = ", E_homo, " ngrid= ",ngrid, " npdos= ",npdos

    X = []
    pDOS = []
    i = 0
    while i<ngrid:
        X.append(E_alp[0] + i * dx)
        
        pdos = []
        j = 0
        while j<npdos:
            pdos.append(0.0)
            j = j + 1
        pDOS.append(pdos)

        i = i + 1

    # so now pDOS[i][j] - contains DOS for j-th projection group (group of atoms perhaps) at i-th grid point
    #============== Now we are ready to make up pDOS ==================
    if(dos_type=="atomic"):

        f1 = open(outfile,"w")

        prefac = 1.0/(var*math.sqrt(2.0*math.pi))
        alp = 0.5/(var*var)

        integ = 0.0 # integral of DOS

        i = 0
        while i<ngrid:

            if emin<=((X[i]-E_homo)) and ((X[i]-E_homo)<=emax):

                j = 0
                line = str(X[i]-E_homo)+"  "
                while j<npdos:
        
                    pDOS[i][j] = 0.0
        
                    #==============================================================
                    ii = 0
                    while ii<Norb:
        
                        wi = 0.0
                        jj=0
                        while jj<Norb:
                            if(At_indx[jj] in specs[j]):
                                wi = wi + MO[ii][jj]**2 # weight of jj-th atomic orbital in ii-th MO
                            jj = jj + 1
        
                        w = prefac*math.exp(-alp*(X[i] - E_alp[ii])**2)  # normalized Gaussian                        
        
                        pDOS[i][j] = pDOS[i][j] + w*wi
        
        
                        ii = ii + 1
                    #==============================================================
                    integ = integ + pDOS[i][j]*dx
        
                    line = line + str(pDOS[i][j])+"  "
                    j = j + 1
                line = line + "\n" 
                f1.write(line)

            i = i + 1

        print "Integral of DOS = ",integ

        f1.close()     

   

   
                           
# Examples of usage:
#gaussian_dos("CD28S38_PM6.LOG","Cd28S38_pm6.txt",[range(0,28),range(28,66)],"atomic",-8.0, 8.0, 0.05,0.025,"unrestricted")
#yaehmop_dos("input_cd28s38_v5.bind.out","Cd28S38_eht_v5.txt",[range(0,28),range(28,66)],"atomic",-10.0, 10.0, 0.05,0.025,"restricted")


gaussian_dos("pop.log","h2oDOS.dat",[range(0,28),range(28,66)],"atomic",-8.0, 8.0, 0.05,0.025,"unrestricted")
