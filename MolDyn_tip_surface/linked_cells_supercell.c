/*
 *  linked_cells_supercell.c
 *  MolDyn_tip_surface
 *
 *  Created by Dimitra Georgakaki on 20/09/2012.
 *  Copyright 2012 Home. All rights reserved.
 */
 
#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void linked_cells_supercell(int mycellsx,int mycellsy, int mycellsz);



//------------------------------------------Build HEAD,LIST and neighbours_supercell for LINKED-CELLS METHOD---------------------------------------------------------------------------------------------------


void linked_cells_supercell(int mycellsx,int mycellsy,int mycellsz)
{
	int index,box,myn; //indexes for arrays
	int ghostsx,ghostsy;
	
	ghostsx = mycellsx+2; //dim of matrix with ghost cells_supercell
	ghostsy = mycellsy+2;
	
	cells_supercell=(int**)malloc((ghostsy)*sizeof(int*)); 
	for (i = 0; i < ghostsy; ++i)
	{cells_supercell[i]=(int*)malloc((ghostsx)*sizeof(int));}
		
	neighbours_supercell=(int**)malloc((26)*sizeof(int*));
	for (i = 0; i < 26; ++i)
	{neighbours_supercell[i]=(int*)malloc((mycellsx*mycellsy*mycellsz)*sizeof(int));} // 26 neighbours_supercell and one central
		
	HEAD_A=(int*)malloc((mycellsx*mycellsy*mycellsz)*sizeof(int));
	HEAD_B=(int*)malloc((mycellsx*mycellsy*mycellsz)*sizeof(int));
	LIST_A=(int*)malloc((npartA)*sizeof(int));
	LIST_B=(int*)malloc((npartB)*sizeof(int));
	
	// Allocate arrays for potential interactions
	partindex_A=(int*)malloc((npartA)*sizeof(int)); // for the interactions in the current cell
	partindex_B=(int*)malloc((npartB)*sizeof(int)); // for the interactions in the current cell
	nextpartindex_B=(int*)malloc((npartB)*sizeof(int)); // for the interactions in the neighbour cells around the current cell

		
//----------------------------------------------Build the neighbours_supercell extended matrix-----------------------------------------------------------------------------
		

		
// We build the cells_supercell matrix 
		
		index=0;
		for (i=ghostsy-2;i>=1;--i){
			for (j=1;j<ghostsx-1;++j){
				index=index+1;
				cells_supercell[i][j]=index;
				
			}
		}
		
		for (j=1;j<=ghostsx-2;++j){
			cells_supercell[0][j]=cells_supercell[ghostsy-2][j];
			cells_supercell[ghostsy-1][j]=cells_supercell[1][j];
			
		}
		
		for (i=0;i<=ghostsy-1;++i){
			cells_supercell[i][ghostsx-1]=cells_supercell[i][1];
			cells_supercell[i][0]=cells_supercell[i][ghostsx-2];
			
		}

		
//We build the neighbours_supercell matrix
		
		for (i=0;i<26;++i){
			for (j=0;j<mycellsx*mycellsy*mycellsz;++j){
				neighbours_supercell[i][j]=0;
			}}
		
		box=0;
		for (i=ghostsy-2;i>=1;--i){
			for (j=1;j<ghostsx-1;++j){
				
				myn=0; //row index
				for (k=-1;k<=1;++k){
					for (l=-1;l<=1;++l){
						if (k!=0 || l!=0)
						{	
							neighbours_supercell[myn][box]=cells_supercell[i+k][j+l];
							myn=myn+1;
						}
					}}
				box=box+1;
				
			}
		}


//We build the neighbours_supercell that are found in parallel planes
	
		for (j=0;j<mycellsx*mycellsy*mycellsz-mycellsx*mycellsy;++j){ 
			for (i=0;i<8;++i){
				
				
				neighbours_supercell[i][j+mycellsx*mycellsy]=neighbours_supercell[i][j]+mycellsx*mycellsy;
				
			}}

	
		
// build the +mycellsx*mycellsy
		
		for (j=0;j<mycellsx*mycellsy*mycellsz;++j){ 
			for (i=0;i<8;++i){
				
				neighbours_supercell[i+8][j]=neighbours_supercell[i][j]+mycellsx*mycellsy;
				
				if (neighbours_supercell[i+8][j]>mycellsx*mycellsy*mycellsz)
					
				{neighbours_supercell[i+8][j]=neighbours_supercell[i][j]-(mycellsz-1)*mycellsx*mycellsy;}
				
				
			}}
		
// build the -mycellsx*mycellsy
		
		for (j=0;j<mycellsx*mycellsy*mycellsz;++j){ 
			for (i=0;i<8;++i){
				
				neighbours_supercell[i+16][j]=neighbours_supercell[i][j] - mycellsx*mycellsy;
				
				if (neighbours_supercell[i+16][j]<=0)
					
				{neighbours_supercell[i+16][j]=neighbours_supercell[i][j]+(mycellsz-1)*mycellsx*mycellsy;}
				
				
			}}
		
//build line 24, the "holes" of the bottom to top layer 
	
		i = 24;
	
		for (j=0;j<mycellsx*mycellsy*mycellsz;++j){
			
			if(j+1 + mycellsx*mycellsy > mycellsx*mycellsy*mycellsz) 
				neighbours_supercell[i][j] = j+1 - (mycellsz - 1)*mycellsx*mycellsy;
			else
				neighbours_supercell[i][j] = j+1 + mycellsx*mycellsy;
		}
	
//build line 25, the "holes" of the top to bottom layer 

		i = 25;
		for (j=0;j<mycellsx*mycellsy*mycellsz;++j){
			if(j+1 - mycellsx*mycellsy <= 0) 
				neighbours_supercell[i][j] = j+1 + (mycellsz - 1)*mycellsx*mycellsy;
			else
				neighbours_supercell[i][j] = j+1 - mycellsx*mycellsy;
		}
		
// compatibility with c
		
		for (i=0;i<26;++i){
			for (j=0;j<mycellsx*mycellsy*mycellsz;++j){
				
				neighbours_supercell[i][j]=neighbours_supercell[i][j]-1;
				
			}}

		


}

