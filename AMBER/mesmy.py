#!/usr/bin/env python3

from argparse import ArgumentParser
import mdtraj as mdt
__version__ = "1.0.0"


def mesmy_cli():
    parser = ArgumentParser(
        description="Create a script for a multi-step Amber MD"
        " relaxation/equilibration workflow."
    )
    parser.add_argument("-i", "--inpcrd",
                        help="Input Amber CRD file.", required=True)
    parser.add_argument("-p", "--prmtop",
                        help="Input Amber PRMTOP file.", required=True)
    parser.add_argument("--version", action="version", version=__version__)

    args = parser.parse_args()

    t = mdt.load(args.inpcrd, top=args.prmtop)
    t_solute = t.atom_slice(
        t.topology.select('not water'))
    nres = t_solute.n_residues
    script = f"""#!/bin/bash

    # An equilibration workflow.
    # Designed for "standard" protein or protein/ligand systems
    # in explicit solvent. Not optimised for membrane protein systems.
    #
    # Developed from scripts provided by the Hughes lab at the University
    # of Montana.
    #

    # You may wish to modify some of the parameters below.
    prmtop_file="{args.prmtop}"
    inpcrd_file="{args.inpcrd}"
    solute={nres} # number of residues in solute
    T="310.0" # target temperature in K
    PMEMD="pmemd.cuda" # name of your MD executable (e.g. may be "pmemd.MPI")

    ### DO NOT MODIFY BELOW THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING ###
    # 1K step energy minimization, strong restraints on heavy atoms, no shake
    cat > step1.in <<EOF
    Min with strong restraints on heavy atoms, no shake
    &cntrl
    imin = 1, ncyc = 50, maxcyc = 1000,
    ntpr = 50,
    ntr = 1, restraintmask = ':1-$solute & !@H=', restraint_wt = 5.0,
    &end
    EOF

    # NVT MD with strong restraints on heavy atoms, shake, dt=.001, 15 ps
    cat > step2.in <<EOF
    NVT MD with strong restraints on heavy atoms, shake dt 0.001 15 ps
    &cntrl
    imin = 0, nstlim = 30000, dt=0.001,
    ntpr = 50, ntwr = 500,
    iwrap = 1,
    ntc = 2, ntf = 2,
    ntt = 3, gamma_ln = 5, temp0 = $T, tempi = $T,
    ntr = 1, restraintmask = ':1-$solute & !@H=', restraint_wt = 5.0,
    &end
    EOF

    # Energy minimization with relaxed restraints on heavy atoms, no shake
    cat > step3.in <<EOF
    Min with relaxed restraints on heavy atoms
    &cntrl
    imin = 1, ncyc = 50, maxcyc = 1000,
    ntpr = 50,
    ntr = 1, restraintmask = ':1-$solute & !@H=', restraint_wt = 2.0,
    &end
    EOF

    # Energy minimization with minimal restraints on heavy atoms, no shake
    cat > step4.in <<EOF
    Min with minimal restraints on heavy atoms
    &cntrl
    imin = 1, ncyc = 50, maxcyc = 1000,
    ntpr = 50,
    ntr = 1, restraintmask = ':1-$solute & !@H=', restraint_wt = 0.1,
    &end
    EOF

    # Energy minimization with no restraints, no shake
    cat > step5.in <<EOF
    Min with no restraints
    &cntrl
    imin = 1, ncyc = 50, maxcyc = 1000,
    ntpr = 50,
    &end
    EOF

    # NPT MD with shake and low restraints on heavy atoms, 20 ps dt=.002
    cat > step6.in <<EOF
    NPT MD with shake and low restraints on heavy atoms, 20 ps dt=.002
    &cntrl
    imin = 0, nstlim = 10000, dt = 0.002,
    ntwx = 1000, ntpr = 50, ntwr = 500,
    iwrap = 1,
    ntc = 2, ntf = 2, ntb = 2,
    ntt = 3, gamma_ln = 5, temp0 = $T, tempi = $T,
    ntp = 1, taup = 1.0,
    ntr = 1, restraintmask = ':1-$solute & !@H=', restraint_wt = 1.0,
    &end
    EOF

    # NPT MD with shake and minimal restraints on heavy atoms
    cat > step7.in <<EOF
    NPT MD with minimal restraints on heavy atoms, 20 ps dt=.002
    &cntrl
    imin = 0, nstlim = 10000, dt=0.002,
    ntx = 5, irest = 1,
    ntwx = 1000, ntpr = 50, ntwr = 500,
    iwrap = 1,
    ntc = 2, ntf = 2, ntb = 2,
    ntt = 3, gamma_ln = 5, temp0 = $T, tempi = $T,
    ntp = 1, taup = 1.0,
    ntr = 1, restraintmask = ':1-$solute & !@H=', restraint_wt = 0.5,
    &end
    EOF

    # NPT MD, shake and minimal restraints on backbone atoms, dt=0.002, 20 ps
    cat > step8.in <<EOF
    NPT MD with minimal restraints on backbone atoms, 20 ps dt=.002
    &cntrl
    imin = 0, nstlim = 10000, dt=0.002,
    ntx = 5, irest = 1,
    ntwx = 1000, ntpr = 50, ntwr = 500,
    iwrap = 1,
    ntc = 2, ntf = 2, ntb = 2,
    ntt = 3, gamma_ln = 5, temp0 = $T, tempi = $T,
    ntp = 1, taup = 1.0,
    ntr = 1, restraintmask = ":1-$solute@H,N,CA,HA,C,O", restraint_wt = 0.5,
    &end
    EOF

    # NPT MD with shake and no restraints, dt=0.002, 200 ps
    cat > step9.in <<EOF
    NPT MD no restraints, dt=0.002, 200 ps
    &cntrl
    imin = 0, nstlim = 10000, dt=0.002,
    ntx = 5, irest = 1,
    ntwx = 1000, ntpr = 50, ntwr = 500,
    iwrap = 1,
    ntc = 2, ntf = 2, ntb = 2,
    ntt = 3, gamma_ln = 5, temp0 = $T, tempi = $T,
    ntp = 1, taup = 1.0,
    &end
    EOF

    START="`date +%s.%N`"

    # Minimization Phase - reference coords are updated each run
    for RUN in step1 step2 step3 step4 step5 ; do
    echo "------------------------"
    echo "Minimization phase: $RUN"
    echo "------------------------"
    if [[ ! -f $RUN.out ]]; then
        echo "File -- $RUN.out -- does not exists. Running job..."
        runme=1
    else
        grep -q 'Total wall time' $RUN.out
        retval=$?
        if [ $retval -ne 0 ]; then
        echo "$RUN did not complete last time. Trying again..."
        runme=1
        else
            echo "$RUN has already completed.  Checking the next step."
            runme=0
        fi
    fi

    if [ $runme -eq 1 ]; then
        $PMEMD -O -i $RUN.in -p $prmtop_file -c $inpcrd_file \
        -ref $inpcrd_file -o $RUN.out -x $RUN.nc -r $RUN.ncrst \
            -inf $RUN.mdinfo
        grep -q 'Total wall time' $RUN.out
        retval=$?
        if [ $retval -ne 0 ]; then
            echo "Error - job failed at this step."
        exit $retval
        fi
    fi

    echo ""
    inpcrd_file="$RUN.ncrst"
    done

    # Equilibration phase: reference coords are last coords from minimize phase
    REF=$inpcrd_file
    for RUN in step6 step7 step8 step9 ; do

    echo "------------------------"
    echo "Equilibration phase: $RUN"
    echo "------------------------"
    if [[ ! -f $RUN.out ]]; then
        echo "File -- $RUN.out -- does not exists. Running job..."
        runme=1
    else
        grep -q 'Total wall time' $RUN.out
        retval=$?
        if [ $retval -ne 0 ]; then
        echo "$RUN did not complete last time. Trying again..."
        runme=1
        else
            echo "$RUN has already completed.  Checking the next step."
            runme=0
        fi
    fi

    if [ $runme -eq 1 ]; then
        $PMEMD -O -i $RUN.in -p $prmtop_file -c $inpcrd_file\
        -ref $REF -o $RUN.out -x $RUN.nc -r $RUN.ncrst -inf $RUN.mdinfo
        grep -q 'Total wall time' $RUN.out
        retval=$?
        if [ $retval -ne 0 ]; then
            echo "Error - job failed at this step."
        exit $retval
        fi
    fi

    echo ""
    inpcrd_file="$RUN.ncrst"
    done

    # Reset the time in the restart file to zero:
    sed -i 's/0.3000000E+02/0.0000000E+00/g' step9.ncrst

    STOP="`date +%s.%N`"
    # The divide by 1 in the line below is to overcome a bug in bc...
    TIMING=`echo "scale=1; ($STOP - $START) / 1;" | bc`
    echo "Total run time: $TIMING seconds."
    echo ""

    exit 0
    """
    print(script)


if __name__ == "__main__":
    mesmy_cli()
