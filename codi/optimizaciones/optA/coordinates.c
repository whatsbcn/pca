/*
This file is part of ftdock, a program for rigid-body protein-protein docking 
Copyright (C) 1997-2000 Gidon Moont

Biomolecular Modelling Laboratory
Imperial Cancer Research Fund
44 Lincoln's Inn Fields
London WC2A 3PX

+44 (0)20 7269 3348
http://www.bmm.icnet.uk/

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/
#include <emmintrin.h>
#include "structures.h"

/************/

int gord( float position , float grid_span , int grid_size ) {

  int ordinate ;

  float one_span = grid_span / (float)grid_size ;

  ordinate = (int)( position / one_span ) + ( grid_size / 2 ) ;

  if( position < 0 ) ordinate -= 1 ;

  return ordinate ;

}

/************/

float pythagoras( float x1 , float y1 , float z1 , float x2 , float y2 , float z2 ) {

	float vector[4];

	//load a registre vectorial	
	__m128 v1 = _mm_set_ps(x1, y1, z1, 0);
	__m128 v2 = _mm_set_ps(x2, y2, z2, 0);

	//operacions
	v1 = _mm_sub_ps(v1, v2);
	v1 = _mm_mul_ps(v1, v1);

	//store de registre vectorial
	_mm_storeu_ps(vector, v1);

	return sqrt(vector[1]+vector[2]+vector[3]);
}
