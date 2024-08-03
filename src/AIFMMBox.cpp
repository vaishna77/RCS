//
//  FMM2DTree.hpp
//
//
//  Created by Vaishnavi Gujjula on 1/4/21.
//
//
#include "AIFMMBox.hpp"

FMM2DBox::FMM2DBox () {
	boxNumber		=	-1;
	parentNumber	=	-1;
	for (int l=0; l<4; ++l) {
		childrenNumbers[l]	=	-1;
	}
	Eliminated = false;
	ILActive = true;
	maxPivot_L2P = 0.0;
	maxPivot_P2M = 0.0;
}
