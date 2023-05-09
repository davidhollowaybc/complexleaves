/*
 *
 *  This file is part of the Virtual Leaf.
 *
 *  The Virtual Leaf is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  The Virtual Leaf is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with the Virtual Leaf.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright 2010 Roeland Merks.
 *
 */

#include <QObject>
#include <QtGui>

#include "simplugin.h"

#include "parameter.h"

#include "wallbase.h"
#include "cellbase.h"
#include "MerksMixed.h"

#include "Pi.h"
#include "random.h"

#include "flux_function.h"

static const std::string _module_id("$Id$");

QString MerksMixed::ModelID(void) {
  // specify the name of your model here
  return QString( "MerksMixed" );
}

// return the number of chemicals your model uses
int MerksMixed::NChem(void) { return 3; }

// To be executed after cell division
void MerksMixed::OnDivide(ParentInfo *parent_info, CellBase *daughter1, CellBase *daughter2) {
  // rules to be executed after cell division go here
  // (e.g., cell differentiation rules)
}

void MerksMixed::SetCellColor(CellBase *c, QColor *color) { 
  // add cell coloring rules here

	double red=c->Chemical(1)/(1.+c->Chemical(1));
	double green=c->Chemical(0)/(1.+c->Chemical(0));
	double blue=0;
//	double blue=c->Chemical(3)/(1.+c->Chemical(3));
	color->setRgbF(red,green,blue);

}

void MerksMixed::CellHouseKeeping(CellBase *c) {
	// add cell behavioral rules here

	  if (c->AtBoundaryP()){

		if(c->Chemical(2) < par->rho0 && c->Chemical(2) > par->rho1){
			c->EnlargeTargetArea(par->mu);
		}
		else if(c->Chemical(2) < par->c0){
			c->EnlargeTargetArea(par->nu);
		}
		else{
		    c->EnlargeTargetArea(par->gamma);
		}

   // 10/23/20 implements aux-threshold for aux-dep growth   if(c->Chemical(0) > par->vessel_inh_level && c->Chemical(2) > 0) c->EnlargeTargetArea(par->cell_expansion_rate*c->Chemical(0)/(par->vessel_expansion_rate + c->Chemical(0)));
		// following: implements Hill kinetics aux-dep growth, no threshold, 11/4/21
		if(c->Chemical(2) > 0) c->EnlargeTargetArea(par->cell_expansion_rate*(pow(c->Chemical(0),par->vessel_inh_level)/(par->vessel_expansion_rate + pow(c->Chemical(0),par->vessel_inh_level))));
	  
	  }
	
	if (c->Area() > 2 * c->BaseArea()) {
		c->Divide();
	} 


}

void MerksMixed::CelltoCellTransport(Wall *w, double *dchem_c1, double *dchem_c2) {
  // add biochemical transport rules here
  
	/* clause from 8/17 to 2/18 to allow limiting patt. form to margin. Removed 2/23/18
		if((w->C2()->AtBoundaryP() && w->C1()->AtBoundaryP()) || (w->C1()->Chemical(0) > par->i1) || (w->C2()->Chemical(0) > par->i1) ) 
			*/
	// fragment of below if condition && w->C1()->CellType() == 0 && w->C2()->CellType() == 0)
	
	// following condition to prevent flow outside of leaf:
	if ((w->C1()->Index() >= 0) && (w->C2()->Index() >= 0) ) {
		  
	    double T = par->transport;
		if ((w->C2()->AtBoundaryP() && !w->C1()->AtBoundaryP()) || (w->C1()->AtBoundaryP() && !w->C2()->AtBoundaryP())) T /= par->c;
	
		// Passive fluxes (Fick's law)

		double phi = w->Length() * par->D[0] * ( w->C2()->Chemical(0) - w->C1()->Chemical(0) );
		dchem_c1[0]+=phi;
		dchem_c2[0]-=phi;


	// directed transport, has saturation of Merks07 eq 1
	// efflux from cell 1 to cell 2
   
 	   double trans12 = ( T * w->Transporters1(1) * 
						  w->C1()->Chemical(0) / (par->ka + w->C1()->Chemical(0)) );
	
    // efflux from cell 2 to cell 1
	    double trans21 = ( T * w->Transporters2(1) * 
						  w->C2()->Chemical(0) / (par->ka + w->C2()->Chemical(0)) );
    
       dchem_c1[0] += trans21 - trans12;
       dchem_c2[0] += trans12 - trans21;
  }
	

	// Influx at leaf "AuxinSource"
	// (as specified in initial condition). zeroth order production of auxin
	 
	if (w->AuxinSource()) { // test if wall is auxin source
	double wall_source = par->leaf_tip_source * w->Length();
	dchem_c1[0] += wall_source;
	dchem_c2[0] += wall_source;
	}  

   // Sink at leaf "AuxinSink"
	// (as specified in initial condition). Changing to 1st order decay of auxin, 8/8/17
	if (w->AuxinSink()) { // test if wall is auxin sink
    //	double wall_sink = par->sam_auxin_breakdown;
	// * w->Length();
	//if(w->C1()->Chemical(0) > wall_sink) dchem_c1[0] -= wall_sink; Old zeroth order removal, pre 8/8/17
	//if(w->C2()->Chemical(0) > wall_sink) dchem_c2[0] -= wall_sink;

	// 10/30/17: changed to biological 'sink', a cell into which auxin flows, via UTG. Therefore added zeroth order production, 
		// as well as first-order decay, to keep steady state. Not keeping 'length' as part of source term. 
		  dchem_c1[0] += par->leaf_tip_source - w->C1()->Chemical(0) * par->sam_auxin_breakdown;
		  dchem_c2[0] += par->leaf_tip_source - w->C2()->Chemical(0) * par->sam_auxin_breakdown;
		  dchem_c1[1] -= par->pin_breakdown_internal * w->C1()->Chemical(1);
		  dchem_c2[1] -= par->pin_breakdown_internal * w->C2()->Chemical(1);
	}  


}
void MerksMixed::WallDynamics(Wall *w, double *dw1, double *dw2) {
  // add biochemical networks for reactions occuring at walls here

	dw1[0] = 0.; dw2[0] = 0.; // chemical 0 unused in walls
	dw1[2] = 0.; dw2[2] = 0.; // chemical 2 unused in walls

 	
   dw1[1] = PINflux(w->C1(),w->C2(),w);
   dw2[1] = PINflux(w->C2(),w->C1(),w); 
  
}
void MerksMixed::CellDynamics(CellBase *c, double *dchem) { 
  // add biochemical networks for intracellular reactions here

// 1st order auxin production in boundary cells    	
	if(c->AtBoundaryP()){
	   if(c->Chemical(2) > par->f){
		 dchem[0] += par->aux1prod * c->Chemical(2) - par->aux1decay * c->Chemical(0);
		// dchem[1] += par->pin_prod * c->Chemical(0) - par->pin_breakdown * c->Chemical(1) - SumFluxFromWalls( c, MerksMixed::PINflux );
	   }

	   else{
		 dchem[0] += (- par->aux1decay * c->Chemical(0)); 
		 dchem[1] += 0.;
	   }

	   if(c->Chemical(2) < par->eps)
	   {
	      dchem[2] += par->aux_cons;
	   }
	}

	else{
		dchem[0] += (par->aux1prodmeso - par->aux1decay * c->Chemical(0));
		
	}

	dchem[1] += par->pin_prod * c->Chemical(0) - par->pin_breakdown * c->Chemical(1) - SumFluxFromWalls( c, MerksMixed::PINflux );
		 
}

 double MerksMixed::PINflux(CellBase *this_cell,
	CellBase *adjacent_cell, Wall *w) { 

   
	double wtf_add = 0.;
	double utg_add = 0.;
	double pin_net = 0.;

  	

 /*  double pin_atwall; // pick the correct side of the Wall
	   if (w->C1() == this_cell) pin_atwall = w->Transporters1(1);
	   else pin_atwall=w->Transporters2(1); */

	double pin_atwall_this, pin_atwall_adj; // pick the correct side of the Wall
	if (w->C1() == this_cell){
		pin_atwall_this = w->Transporters1(1);
		pin_atwall_adj = w->Transporters2(1);
	}
	else {
		pin_atwall_this=w->Transporters2(1);
		pin_atwall_adj=w->Transporters1(1);
	}

	 // calculate PIN translocation rate from cell to membrane; WTF (R-L and P, 05), but retaining saturated T flow of Merks07 eq1. 
	  // This calculates phi_tot, total flux, in order to have a current value for the wall allocation step. 
	// fragment of below condition, && w->C1()->CellType() == 0 && w->C2()->CellType() == 0)
	
	if ((w->C1()->Index() >= 0) && (w->C2()->Index() >= 0))  {
		
		double T = par->transport;
		if ((w->C2()->AtBoundaryP() && !w->C1()->AtBoundaryP()) || (w->C1()->AtBoundaryP() && !w->C2()->AtBoundaryP())) T /= par->c;

		double trans_thisadj = (T * pin_atwall_this * this_cell->Chemical(0) / (par->ka + this_cell->Chemical(0)));

		double trans_adjthis = (T * pin_atwall_adj * adjacent_cell->Chemical(0) / (par->ka + adjacent_cell->Chemical(0)));

		double phi_tot = (w->Length()*par->D[0] * (this_cell->Chemical(0) - adjacent_cell->Chemical(0)))
			+ trans_thisadj - trans_adjthis;


		if (phi_tot > 0){
//			if ((this_cell->AtBoundaryP() && this_cell->Chemical(0) > par->f) || (!this_cell->AtBoundaryP())) { USES OLD par->f, pre 11/23/17
			// next line added 2-8-18. par->f threshold for wtf add rm-ed 5/10/18, but it's been set arbitrarily high prior to this. 
				wtf_add = (this_cell->Chemical(1) / (par->kap + this_cell->Chemical(1))) * (par->e * pow(phi_tot, 2) + par->d * phi_tot);
		}
	}
	  	
		// version from R-L + P 05, though they had a hard stop of Pi<4 (no ij in their model)
		//	pin_flux = par->k1 * pow(phi_tot, 2) + par->d - par->k2 * pin_atwall1; 
	    
	// UTG allocation, added back in, 10/16/17. AtBoundary implements only for boundary cells. 
	
      double receptor_level = adjacent_cell->Chemical(0) * par->r / (par->kr + adjacent_cell->Chemical(0));  
	      
	     utg_add = par->k1 * this_cell->Chemical(1) * receptor_level / ( par->km + this_cell->Chemical(1) );

	// total allocation

		   if(this_cell->AtBoundaryP() && adjacent_cell->AtBoundaryP()){ 
			   if(this_cell->CellType() == adjacent_cell->CellType()){
					pin_net = utg_add + wtf_add - par->k2 * pin_atwall_this;
		       }
			   else{
				   pin_net = 0.;
			   }
		   }


		   else {
			   	   pin_net = utg_add + wtf_add - par->k2 * pin_atwall_this;
		   }
	
	// rm-ed 12/2/19, nu had been 2000, i.e. Pij 'ceiling' not used	 if (pin_atwall_this > par->nu) pin_net = 0.0;
		  
  return pin_net;  


	
}


Q_EXPORT_PLUGIN2(MerksMixed, MerksMixed)
