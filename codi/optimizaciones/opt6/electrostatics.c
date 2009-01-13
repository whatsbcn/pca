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

#include "xmmintrin.h"
#include "structures.h"

void assign_charges( struct Structure This_Structure ) {

/************/

  /* Counters */

  int	residue , atom ;

/************/

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      This_Structure.Residue[residue].Atom[atom].charge = 0.0 ;

      /* peptide backbone */

      if( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " N  " ) == 0 ) {
        if( strcmp( This_Structure.Residue[residue].res_name , "PRO" ) == 0 ) {
          This_Structure.Residue[residue].Atom[atom].charge = -0.10 ;
        } else {
          This_Structure.Residue[residue].Atom[atom].charge =  0.55 ;
          if( residue == 1 ) This_Structure.Residue[residue].Atom[atom].charge = 1.00 ;
        }
      }

      if( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " O  " ) == 0 ) {
        This_Structure.Residue[residue].Atom[atom].charge = -0.55 ;
        if( residue == This_Structure.length  ) This_Structure.Residue[residue].Atom[atom].charge = -1.00 ;
      }

      /* charged residues */

      if( ( strcmp( This_Structure.Residue[residue].res_name , "ARG" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " NH" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge =  0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "ASP" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " OD" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge = -0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "GLU" ) == 0 ) && ( strncmp( This_Structure.Residue[residue].Atom[atom].atom_name , " OE" , 3 ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge = -0.50 ;
      if( ( strcmp( This_Structure.Residue[residue].res_name , "LYS" ) == 0 ) && ( strcmp( This_Structure.Residue[residue].Atom[atom].atom_name , " NZ " ) == 0 ) ) This_Structure.Residue[residue].Atom[atom].charge =  1.00 ;

    }
  }

/************/

}



/************************/



void electric_field( struct Structure This_Structure , float grid_span , int grid_size , fftw_real *grid ) {

/************/

  /* Counters */

  int	residue , atom ;

  /* Co-ordinates */

  int	x , y , z , i , j , num_unrolled_iters;
  float		x_centre , y_centre , z_centre ;

  /* Variables */

  float		distance[4];
  float		phi[4] __attribute__((aligned(16))) , epsilon[4] , coefficient[4] __attribute__((aligned(16)));

  __m128 	phiVect , coefficientVect , atomVect;

  char print_buffer[grid_size];

  int indexCharge , indexCoord;

  float charge[4400] __attribute__((aligned(16)));
  float coord[4400*3];

/************/

  i = 0;
  j = 0;
  while( j < grid_size*grid_size*grid_size) {
    for( i = j ; i < j+grid_size ; i ++ ) {
      grid[i] = (fftw_real)0;
    }
    j += grid_size + 2;
  }

  indexCoord = 0;
  indexCharge = 0;
  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {
      if( This_Structure.Residue[residue].Atom[atom].charge != 0 ) {
        charge[indexCharge] = This_Structure.Residue[residue].Atom[atom].charge;
        coord[indexCoord] = This_Structure.Residue[residue].Atom[atom].coord[1];
        coord[indexCoord+1] = This_Structure.Residue[residue].Atom[atom].coord[2];
        coord[indexCoord+2] = This_Structure.Residue[residue].Atom[atom].coord[3];
        indexCharge++;
        indexCoord+=3;
      }
    }
  }

/************/

  setvbuf( stdout , (char *)NULL , _IONBF , 0 ) ;

  printf( "  electric field calculations ( one dot / grid sheet ) " ) ;

  i = 0;
  for( x = 0 ; x < grid_size ; x ++ ) {

    print_buffer[x] = '.';
//    printf(".");

    x_centre  = gcentre( x , grid_span , grid_size ) ;

    for( y = 0 ; y < grid_size ; y ++ ) {

      y_centre  = gcentre( y , grid_span , grid_size ) ;

      for( z = 0 ; z < grid_size ; z ++ ) {

        z_centre  = gcentre( z , grid_span , grid_size ) ;

        phiVect = _mm_set1_ps(0.0);

        indexCoord = 0;
        num_unrolled_iters = indexCharge - (indexCharge % 4);

        for( atom = 0 ; atom <= num_unrolled_iters ; atom += 4 ) {

          distance[0] = pythagoras( coord[indexCoord], coord[indexCoord+1], coord[indexCoord+2] , x_centre , y_centre , z_centre ) ;
          distance[1] = pythagoras( coord[indexCoord+3], coord[indexCoord+4], coord[indexCoord+5] , x_centre , y_centre , z_centre ) ;
          distance[2] = pythagoras( coord[indexCoord+6], coord[indexCoord+7], coord[indexCoord+8] , x_centre , y_centre , z_centre ) ;
          distance[3] = pythagoras( coord[indexCoord+9], coord[indexCoord+10], coord[indexCoord+11] , x_centre , y_centre , z_centre ) ;
 
          indexCoord += 12;

          atomVect = _mm_load_ps(&(charge[atom]));

          if( distance[0] < 2.0 ) distance[0] = 2.0 ;
          if( distance[1] < 2.0 ) distance[1] = 2.0 ;
          if( distance[2] < 2.0 ) distance[2] = 2.0 ;
          if( distance[3] < 2.0 ) distance[3] = 2.0 ;

          if (distance[0] >= 8.0)
            coefficient[0] = distance[0] * 80.0;
          else if (distance[0] <= 6.0)
            coefficient[0] = distance[0] * 4.0;
          else
            coefficient[0] = (38 * distance[0] - 224) * distance[0];

          if (distance[1] >= 8.0)
            coefficient[1] = distance[1] * 80.0;
          else if (distance[1] <= 6.0)
            coefficient[1] = distance[1] * 4.0;
          else
            coefficient[1] = (38 * distance[1] - 224) * distance[1];

          if (distance[2] >= 8.0)
            coefficient[2] = distance[2] * 80.0;
          else if (distance[2] <= 6.0)
            coefficient[2] = distance[2] * 4.0;
          else
            coefficient[2] = (38 * distance[2] - 224) * distance[2];

          if (distance[3] >= 8.0)
            coefficient[3] = distance[3] * 80.0;
          else if (distance[3] <= 6.0)
            coefficient[3] = distance[3] * 4.0;
          else
            coefficient[3] = (38 * distance[3] - 224) * distance[3];

          coefficientVect = _mm_load_ps(coefficient);

          atomVect = _mm_div_ps(atomVect,coefficientVect);

          phiVect = _mm_add_ps(phiVect,atomVect); 
        }

        _mm_store_ps(phi,phiVect);

        for( atom ; atom <= indexCharge ; atom ++ ) {

          distance[0] = pythagoras( coord[indexCoord], coord[indexCoord+1], coord[indexCoord+2] , x_centre , y_centre , z_centre ) ;

          indexCoord += 3;

          if( distance[0] < 2.0 ) distance[0] = 2.0 ;

          if (distance[0] >= 8.0)
            coefficient[0] = distance[0] * 80.0;
          else if (distance[0] <= 6.0)
            coefficient[0] = distance[0] * 4.0;
          else
            coefficient[0] = (38 * distance[0] - 224) * distance[0];

          phi[0] += charge[atom] / coefficient[0] ;
        }

        phi[0] = phi[0] + phi[1] + phi[2] + phi[3];

        grid[i] = (fftw_real)(phi[0]);
        i++;
      }
      i+=2;
    }
  }

  printf("%s\n",print_buffer);
//  printf("\n");

/************/

  return ;

}



/************************/



void electric_point_charge( struct Structure This_Structure , float grid_span , int grid_size , fftw_real *grid ) {

/************/

  /* Counters */

  int	residue , atom ;

  /* Co-ordinates */

  int	x , y , z , i , j;
  int	x_low , x_high , y_low , y_high , z_low , z_high ;

  float		a , b , c ;
  float		x_corner , y_corner , z_corner ;
  float		w ;

  /* Variables */

  float		one_span ;

/************/

  i = 0;
  j = 0;
  while( j < grid_size*grid_size*grid_size) {
    for( i = j ; i < j+grid_size ; i ++ ) {
      grid[i] = (fftw_real)0;
    }
    j += grid_size + 2;
  }

/************/

  one_span = grid_span / (float)grid_size ;

  for( residue = 1 ; residue <= This_Structure.length ; residue ++ ) {
    for( atom = 1 ; atom <= This_Structure.Residue[residue].size ; atom ++ ) {

      if( This_Structure.Residue[residue].Atom[atom].charge != 0 ) {

        x_low = gord( This_Structure.Residue[residue].Atom[atom].coord[1] - ( one_span / 2 ) , grid_span , grid_size ) ;
        y_low = gord( This_Structure.Residue[residue].Atom[atom].coord[2] - ( one_span / 2 ) , grid_span , grid_size ) ;
        z_low = gord( This_Structure.Residue[residue].Atom[atom].coord[3] - ( one_span / 2 ) , grid_span , grid_size ) ;

        x_high = x_low + 1 ;
        y_high = y_low + 1 ;
        z_high = z_low + 1 ;

        a = This_Structure.Residue[residue].Atom[atom].coord[1] - gcentre( x_low , grid_span , grid_size ) - ( one_span / 2 ) ;
        b = This_Structure.Residue[residue].Atom[atom].coord[2] - gcentre( y_low , grid_span , grid_size ) - ( one_span / 2 ) ;
        c = This_Structure.Residue[residue].Atom[atom].coord[3] - gcentre( z_low , grid_span , grid_size ) - ( one_span / 2 ) ;

        for( x = x_low ; x <= x_high  ; x ++ ) {
 
          x_corner = one_span * ( (float)( x - x_high ) + .5 ) ;

          for( y = y_low ; y <= y_high  ; y ++ ) {

            y_corner = one_span * ( (float)( y - y_high ) + .5 ) ;

            for( z = z_low ; z <= z_high  ; z ++ ) {

              z_corner = one_span * ( (float)( z - z_high ) + .5 ) ;

              w = ( ( x_corner + a ) * ( y_corner + b ) * ( z_corner + c ) ) / ( 8.0 * x_corner * y_corner * z_corner ) ;

              grid[gaddress(x,y,z,grid_size)] += (fftw_real)( w * This_Structure.Residue[residue].Atom[atom].charge ) ;

            }
          }
        }

      }

    }
  }

/************/

  return ;

}



/************************/



void electric_field_zero_core( int grid_size , fftw_real *elec_grid , fftw_real *surface_grid , float internal_value ) {

/************/

  /* Co-ordinates */

  int	x , y , z , i , j;

/************/

  i = 0;
  j = 0;
  while( j < grid_size*grid_size*grid_size) {
    for( i = j ; i < j+grid_size ; i ++ ) {
      if( surface_grid[i] == (fftw_real)internal_value ) elec_grid[i] = (fftw_real)0 ;
    }
    j += grid_size + 2;
  }
/************/

  return ;

}
