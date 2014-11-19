#!/bin/bash

VERBOSITY="ERROR"
GLOBALFLAGS="nopred"

# run all extration versions for yield:
../bin/extract -t yield -i /nfs/dust/cms/user/spanns/ResultsPseudo172OldKinReco_2_finebinning -v $VERBOSITY -s -f $GLOBALFLAGS,nonorm -o /nfs/dust/cms/user/spanns/ExtractedPlotsPAS2/yield-fine

../bin/extract -t yield -i /nfs/dust/cms/user/spanns/ResultsPseudo172OldKinReco_2_finebinning -v $VERBOSITY -s -f $GLOBALFLAGS,bgr,nonorm -o /nfs/dust/cms/user/spanns/ExtractedPlotsPAS2/yield-fine-bgr

../bin/extract -t yield -i /nfs/dust/cms/user/spanns/ResultsPseudo172OldKinReco_2_finebinning -v $VERBOSITY -s -f $GLOBALFLAGS,nopred,norm -o /nfs/dust/cms/user/spanns/ExtractedPlotsPAS2/yield-fine-norm

../bin/extract -t yield -i /nfs/dust/cms/user/spanns/ResultsPseudo172OldKinReco_2_finebinning -v $VERBOSITY -s -f $GLOBALFLAGS,norm,bgr -o /nfs/dust/cms/user/spanns/ExtractedPlotsPAS2/yield-fine-norm-bgr
