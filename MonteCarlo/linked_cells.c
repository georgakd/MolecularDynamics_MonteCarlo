/*
 *  linked_cells.c
 *  tersoff calculations
 *
 *  Created by Dimitra Georgakaki on 06/12/2011.
 *  Copyright 2011 Home. All rights reserved.
 */
 

#include "global.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void linked_cells(int cellsx,int cellsy, int cellsz);



//------------------------------------------Build HEAD,LIST and neighbours for LINKED-CELLS METHOD---------------------------------------------------------------------------------------------------


void linked_cells(int cellsx,int cellsy,int cellsz)
{

	totalcells = cellsx*cellsy*cellsz; // cells in 3D	
	ghostsx = cellsx+2; //dim of matrix with ghost cells
	ghostsy = cellsy+2;

		
	cells=(int**)malloc((ghostsy)*sizeof(int*)); 
	for (i = 0; i < ghostsy; ++i)
	{cells[i]=(int*)malloc((ghostsx)*sizeof(int));}
	
	
	neighbours=(int**)malloc((26)*sizeof(int*));
	for (i = 0; i < 26; ++i)
	{neighbours[i]=(int*)malloc((totalcells)*sizeof(int));} // 26 neighbours and one central
	
	
	HEAD=(int*)malloc((totalcells)*sizeof(int));
	LIST=(int*)malloc((npart)*sizeof(int));
	
	// Allocate arrays for potential interactions

		if (mypotential==0){
		cellpartarray = (int*)malloc((npart)*sizeof(int)); //1D array to place the particle indices of the active central cell
	
		neighborcellpartarray = (int*)malloc((npart)*sizeof(int)); //1D array to place the particle indices of the active neighboring cell
		allpartarray = (int*)malloc((npart)*sizeof(int)); //(*NEW! 1D array to place the particle indices of the central cell and ALL the neighbor cells*)
		}
		if (mypotential==1){
		partindex=(int*)malloc((npart)*sizeof(int)); // for the interactions in the current cell
		nextpartindex=(int*)malloc((npart)*sizeof(int)); // for the interactions in the neighbour cells around the current cell
		}

		
//----------------------------------------------Build the neighbours extended matrix-----------------------------------------------------------------------------
		

		
// We build the cells matrix 
		
		myindex=0;
		for (i=ghostsy-2;i>=1;--i){
			for (j=1;j<ghostsx-1;++j){
				myindex=myindex+1;
				cells[i][j]=myindex;
				
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
			for (j=0;j<totalcells;++j){
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
	
		for (j=0;j<totalcells-cellsx*cellsy;++j){ 
			for (i=0;i<8;++i){
				
				
				neighbours[i][j+cellsx*cellsy]=neighbours[i][j]+cellsx*cellsy;
				
			}}
		
		
// build the +cellsx*cellsy
		
		for (j=0;j<totalcells;++j){ 
			for (i=0;i<8;++i){
				
				neighbours[i+8][j]=neighbours[i][j]+cellsx*cellsy;
				
				if (neighbours[i+8][j]>totalcells)
					
				{neighbours[i+8][j]=neighbours[i][j]-(cellsz-1)*cellsx*cellsy;}
				
				
			}}
		
// build the -cellsx*cellsy
		
		for (j=0;j<totalcells;++j){ 
			for (i=0;i<8;++i){
				
				neighbours[i+16][j]=neighbours[i][j] - cellsx*cellsy;
				
				if (neighbours[i+16][j]<=0)
					
				{neighbours[i+16][j]=neighbours[i][j]+(cellsz-1)*cellsx*cellsy;}
				
				
			}}
		
//build line 24, the "holes" of the bottom to top layer 
	
		i = 24;
	
		for (j=0;j<totalcells;++j){
			
			if(j+1 + cellsx*cellsy > totalcells) 
				neighbours[i][j] = j+1 - (cellsz - 1)*cellsx*cellsy;
			else
				neighbours[i][j] = j+1 + cellsx*cellsy;
		}
	
//build line 25, the "holes" of the top to bottom layer 

		i = 25;
		for (j=0;j<totalcells;++j){
			if(j+1 - cellsx*cellsy <= 0) 
				neighbours[i][j] = j+1 + (cellsz - 1)*cellsx*cellsy;
			else
				neighbours[i][j] = j+1 - cellsx*cellsy;
		}
		
// compatibility with c
		
		for (i=0;i<26;++i){
			for (j=0;j<totalcells;++j){
				
				neighbours[i][j]=neighbours[i][j]-1;
				
			}}
		


}

