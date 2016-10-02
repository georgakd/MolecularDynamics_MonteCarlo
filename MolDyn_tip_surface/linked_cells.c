/*
 *  linked_cells.c
 *  MolDyn_tip_surface
 *
 *  Created by Dimitra Georgakaki on 20/09/2012.
 *  Copyright 2012 Home. All rights reserved.
 */
 

#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void linked_cells(int mycellsx,int mycellsy, int mycellsz);



//------------------------------------------Build HEAD,LIST and neighbours for LINKED-CELLS METHOD---------------------------------------------------------------------------------------------------


void linked_cells(int mycellsx,int mycellsy,int mycellsz)
{
	int index,box,myn; //indexes for arrays
	int ghostsx,ghostsy;
		
	ghostsx = mycellsx+2; //dim of matrix with ghost cells
	ghostsy = mycellsy+2;
		
	cells=(int**)malloc((ghostsy)*sizeof(int*)); 
	for (i = 0; i < ghostsy; ++i)
	{cells[i]=(int*)malloc((ghostsx)*sizeof(int));}
		
	neighbours=(int**)malloc((26)*sizeof(int*));
	for (i = 0; i < 26; ++i)
	{neighbours[i]=(int*)malloc((mycellsx*mycellsy*mycellsz)*sizeof(int));} // 26 neighbours and one central
		
	HEAD=(int*)malloc((mycellsx*mycellsy*mycellsz)*sizeof(int));
	LIST=(int*)malloc((npartB)*sizeof(int));
	
	// Allocate arrays for potential interactions

		if (mypotential==0){ //Tersoff
		cellpartarray = (int*)malloc((npartB)*sizeof(int)); //1D array to place the particle indices of the active central cell	
		neighborcellpartarray = (int*)malloc((npartB)*sizeof(int)); //1D array to place the particle indices of the active neighboring cell
		allpartarray = (int*)malloc((npartB)*sizeof(int)); //(*NEW! 1D array to place the particle indices of the central cell and ALL the neighbor cells*)
		}
		if (mypotential==1){ //LJ
		partindex=(int*)malloc((npartB)*sizeof(int)); // for the interactions in the current cell
		nextpartindex=(int*)malloc((npartB)*sizeof(int)); // for the interactions in the neighbour cells around the current cell
		}

		
//----------------------------------------------Build the neighbours extended matrix-----------------------------------------------------------------------//

		
// We build the cells matrix 
		
		index=0;
		for (i=ghostsy-2;i>=1;--i){
			for (j=1;j<ghostsx-1;++j){
				index=index+1;
				cells[i][j]=index;
				
			}
		}
		
		for (j=1;j<=ghostsx-2;++j){
			cells[0][j]=cells[ghostsy-2][j];
			cells[ghostsy-1][j]=cells[1][j];
			
		}
		
		for (i=0;i<=ghostsy-1;++i){
			cells[i][ghostsx-1]=cells[i][1];
			cells[i][0]=cells[i][ghostsx-2];
			
		}

		
//We build the neighbours matrix
		
		for (i=0;i<26;++i){
			for (j=0;j<mycellsx*mycellsy*mycellsz;++j){
				neighbours[i][j]=0;
			}}
		
		box=0;
		for (i=ghostsy-2;i>=1;--i){
			for (j=1;j<ghostsx-1;++j){
				
				myn=0; //row index
				for (k=-1;k<=1;++k){
					for (l=-1;l<=1;++l){
						if (k!=0 || l!=0)
						{
							
							neighbours[myn][box]=cells[i+k][j+l];
							myn=myn+1;
						}
					}}
				box=box+1;
				
			}
		}

//We build the neighbours that are found in parallel planes
	
		for (j=0;j<mycellsx*mycellsy*mycellsz-mycellsx*mycellsy;++j){ 
			for (i=0;i<8;++i){
				
				
				neighbours[i][j+mycellsx*mycellsy]=neighbours[i][j]+mycellsx*mycellsy;
				
			}}
		
		
// build the +mycellsx*mycellsy
		
		for (j=0;j<mycellsx*mycellsy*mycellsz;++j){ 
			for (i=0;i<8;++i){
				
				neighbours[i+8][j]=neighbours[i][j]+mycellsx*mycellsy;
				
				if (neighbours[i+8][j]>mycellsx*mycellsy*mycellsz)
					
				{neighbours[i+8][j]=neighbours[i][j]-(mycellsz-1)*mycellsx*mycellsy;}
				
				
			}}
		
// build the -mycellsx*mycellsy
		
		for (j=0;j<mycellsx*mycellsy*mycellsz;++j){ 
			for (i=0;i<8;++i){
				
				neighbours[i+16][j]=neighbours[i][j] - mycellsx*mycellsy;
				
				if (neighbours[i+16][j]<=0)
					
				{neighbours[i+16][j]=neighbours[i][j]+(mycellsz-1)*mycellsx*mycellsy;}
				
				
			}}
		
//build line 24, the "holes" of the bottom to top layer 
	
		i = 24;
	
		for (j=0;j<mycellsx*mycellsy*mycellsz;++j){
			
			if(j+1 + mycellsx*mycellsy > mycellsx*mycellsy*mycellsz) 
				neighbours[i][j] = j+1 - (mycellsz - 1)*mycellsx*mycellsy;
			else
				neighbours[i][j] = j+1 + mycellsx*mycellsy;
		}
	
//build line 25, the "holes" of the top to bottom layer 

		i = 25;
		for (j=0;j<mycellsx*mycellsy*mycellsz;++j){
			if(j+1 - mycellsx*mycellsy <= 0) 
				neighbours[i][j] = j+1 + (mycellsz - 1)*mycellsx*mycellsy;
			else
				neighbours[i][j] = j+1 - mycellsx*mycellsy;
		}
		
// compatibility with c
		
		for (i=0;i<26;++i){
			for (j=0;j<mycellsx*mycellsy*mycellsz;++j){
				
				neighbours[i][j]=neighbours[i][j]-1;
				
			}}
		


}

