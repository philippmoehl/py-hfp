/****************************************************************************
* EfPiSoft                                                                  *
*                                                                           *
* Consiglio Nazionale delle Ricerche                                        *
* Istituto di Matematica Applicata e Tecnologie Informatiche                *
* Sezione di Genova                                                         *
* IMATI-GE / CNR                                                            *
*                                                                           *
* Authors: Marco Attene                                                     *
*                                                                           *
* Copyright(C) 2006: IMATI-GE / CNR                                         *
*                                                                           *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#ifndef FITTING_PRIMITIVES_H
#define FITTING_PRIMITIVES_H

#ifdef WIN32
#define SOQT_NOT_DLL
#endif

#include "jmesh.h"
#include "clusterGraph.h"
#include <stdlib.h>

// These are constants to be used within the mask 'what_to_fit'
#define HFP_FIT_PLANES    ((unsigned char)1)
#define HFP_FIT_SPHERES   ((unsigned char)2)
#define HFP_FIT_CYLINDERS ((unsigned char)4)

//////////////////////////////////////////////////////////////////////////
//
// Node of the cluster graph
//
//////////////////////////////////////////////////////////////////////////

class curvcNode : public graphNode
{
public:

	static List* collapses;

	int id;                // Unique identifier of the node
	List triangles;        // All the triangles within the cluster
	List areas;            // Triangle areas
	Point sum_ctr;         // Weighted sum of barycenters
	SymMatrix3x3 Cov_v;    // Covariance matrix of cluster vertices
	SymMatrix3x3 Cov_c;    // Covariance matrix of normal variation
	SymMatrix4x4 AtA;      // To find bf-sphere (matrix)
	double Atb[4];         // To find bf-sphere (known term)
	double tot_area;       // Total area of the cluster

	curvcNode(Triangle*, int);	// Constructor
	~curvcNode() { areas.freeNodes(); }

	static void merge(const void* n1, const void* n2);
	static double edgeCostFunction(const void* n1, const void* n2);

	static bool bestFittingCircle(double*, int, double&, double&, double&);
	static bool bestFittingSphere(SymMatrix4x4&, double*, double&, double&, double&, double&);

	static double fittingPlaneCost(const void*, const void*);
	static double fittingSphereCost(const void*, const void*);
	static double fittingCylinderCost(const void*, const void*);

	static unsigned char what_to_fit; // Binary mask defining fitting primitives

	static bool is_fitting_planes() { return (what_to_fit & HFP_FIT_PLANES); }
	static bool is_fitting_spheres() { return (what_to_fit & HFP_FIT_SPHERES); }
	static bool is_fitting_cylinders() { return (what_to_fit & HFP_FIT_CYLINDERS); }
};

class HFP_Action
{
	Triangulation* tin;

public:

	HFP_Action(Triangulation* t) { tin = t; }

	List* fit(unsigned char what_to_fit);
};

#endif // FITTING_PRIMITIVES_H