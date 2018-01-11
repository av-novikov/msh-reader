#include <string>
#include <iostream>
#include <memory>

#include "src/Scene.hpp"
#include "src/models/oil/Oil.hpp"
#include "src/models/oil/OilSolver.hpp"

using namespace std;

oil::Properties* getProps()
{
	oil::Properties* props = new oil::Properties();
	props->timePeriods.push_back(86400.0 * 20.0);

	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	//props->rates.push_back(100.0);
	props->pwf.push_back(200.0 * 100000.0);
	props->ht = 1000.0;
	props->ht_min = 1000.0;
	props->ht_max = 100000.0;

	props->perfIntervals.push_back(make_pair(0, 0));
	props->R_dim = 10.0;
	props->r_w = 10;
	props->r_e = 1000.0;

	oil::Skeleton_Props tmp;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 200.0 * 100000.0;
	tmp.height = 10.0;
	tmp.kx = tmp.ky = 10.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 0.0 * 1.0e-10;
	props->props_sk.push_back(tmp);

	props->props_oil.beta = 1.e-9;
	props->props_oil.dens_stc = 800.0;
	props->props_oil.visc = 1.0;
	props->props_oil.p_ref = tmp.p_init;

	props->alpha = 7200.0;

	return props;
}

int main()
{
	const auto props = getProps();
	Scene<oil::Oil, oil::OilSolver, oil::Properties> scene;
	scene.load(*props, "attempt4.nebr");
	scene.start();

	return 0;
}