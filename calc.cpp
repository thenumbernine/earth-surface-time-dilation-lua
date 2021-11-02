#include <iostream>
#include <cmath>
#include <tuple>
#include <vector>
#include <functional>
#include <fstream>

template<typename T, typename U>
std::vector<U> map(std::vector<T> const & x, std::function<U(T)> f) {
	std::vector<U> y;
	std::transform(x.cbegin(), x.cend(), std::back_inserter(y), f); 
	return y;
}

int main() {
	using real = double;

	real speedOfLight_in_m_per_s = 299792458;
	std::cout << "speed of light, in m/s: " << speedOfLight_in_m_per_s << std::endl;

	real gravitationalConstant_in_m3_per_kg_s2 = 6.67384e-11;
	std::cout << "gravitational constant, in m^3 / (kg s^2): " << gravitationalConstant_in_m3_per_kg_s2 << std::endl;

	real earth_equatorialRadius_in_m = 6378136;
	std::cout << "Earth equitorial radius in m: " << earth_equatorialRadius_in_m << std::endl;

	real earth_inverseFlattening = 298.257223563;
	std::cout << "Earth inverse flattening: " << earth_inverseFlattening << std::endl;

	real earthRotPeriod_in_s = (23 + (56 + 4.09053083288 / 60) / 60) * 60 * 60;
	std::cout << "Earth rotation period in s: " << earthRotPeriod_in_s << std::endl;

	real earthMass_in_kg = 5.9736e+24;
	std::cout << "Earth mass in kg: " << earthMass_in_kg << std::endl;

	auto inverseFlattening = earth_inverseFlattening;
	auto equatorialRadius_in_m = earth_equatorialRadius_in_m;
	auto rotationPeriod_in_s = earthRotPeriod_in_s;
	auto mass_in_kg = earthMass_in_kg;

	auto flattening = 1 / inverseFlattening;
	auto eccentricitySquared = flattening * (2 - flattening);
	auto angularVelocity_in_1_per_s = (2 * M_PI) / rotationPeriod_in_s;
	std::cout << "angularVelocity_in_1_per_s " << angularVelocity_in_1_per_s << std::endl;

	auto schwarzschildRadius_in_m = 2 * mass_in_kg * gravitationalConstant_in_m3_per_kg_s2 / (speedOfLight_in_m_per_s * speedOfLight_in_m_per_s);

	auto degToRad = [](auto deg) { return deg * M_PI / 180; };

	// lat lon in degrees
	auto geodeticPos = [&](real lat, real lon, real height) {
		auto phi = degToRad(lat);
		auto lambda = degToRad(lon);
		auto cosPhi = cos(phi);
		auto sinPhi = sin(phi);
		auto eccentricitySquared = (2 * inverseFlattening - 1) / (inverseFlattening * inverseFlattening);
		auto sinPhiSquared = sinPhi * sinPhi;
		auto N = equatorialRadius_in_m / sqrt(1 - eccentricitySquared * sinPhiSquared);
		auto NPlusH = N + height;
		auto x = NPlusH * cosPhi * cos(lambda);
		auto y = NPlusH * cosPhi * sin(lambda);
		auto z = (N * (1 - eccentricitySquared) + height) * sinPhi;
		return std::make_tuple(x, y, z);	// in meters
	};

	auto geodeticDist = [&](real lat, real lon, real height) {
		auto [x,y,z] = geodeticPos(lat, lon, height);
		return sqrt(x*x + y*y + z*z);
	};

	auto geodeticVel = [&](real lat, real lon, real height) {
		auto phi = degToRad(lat);
		auto lambda = degToRad(lon);
		auto cosPhi = cos(phi);
		auto sinPhi = sin(phi);
		auto eccentricitySquared = (2 * inverseFlattening - 1) / (inverseFlattening * inverseFlattening);
		auto sinPhiSquared = sinPhi * sinPhi;
		auto N = equatorialRadius_in_m / sqrt(1 - eccentricitySquared * sinPhiSquared);
		auto NPlusH = N + height;
		auto x = -angularVelocity_in_1_per_s * NPlusH * cosPhi * sin(lambda);
		auto y = angularVelocity_in_1_per_s * NPlusH * cosPhi * cos(lambda);
		auto z = 0;
		return std::make_tuple(x, y, z);
	};

	auto geodeticSpeed = [&](real lat, real lon, real height) {
		auto [x, y, z] = geodeticVel(lat, lon, height);
		return sqrt(x*x + y*y + z*z);
	};

	constexpr int n = 180;
	std::vector<real> x;
	for (int i = 0; i < n; ++i) {
		x.push_back( ((real)i - 1) / ((real)n - 1) );
	}
	auto lats = map<real, real>(x, [](real x){ return x * 180 - 90; });
	auto dists = map<real, real>(lats, [&](real lat){ return geodeticDist(lat, 0, 0); });
	auto vels = map<real, real>(lats, [&](real lat){ return geodeticSpeed(lat, 0, 0); });
	auto lorentzBetas = map<real, real>(vels, [&](real v) { return v / speedOfLight_in_m_per_s; });
	auto lorentzGammas = map<real, real>(lorentzBetas, [&](real beta){ return sqrt(1 - beta * beta); });
	/*
	Schwarzschild metric: ds^2 = -(1 - R/r) dt^2 + 1/(1 - R/r) dr^2 + r^2 dθ^2 + r^2 sin(θ)^2 dφ^2
	for an object not moving radially, not moving latitudally ...
	ds^2 = -(1 - R/r) dt^2 + r^2 sin(θ)^2 dφ^2 
	... and dλ/ds = omega, right? or is it dλ/dt? 
	(ds/dt)^2 = -(1 - R/r) + r^2 sin(θ)^2 (dφ/dt)^2 
	... in geodetic coordinates ...
	(ds/dt)^2 = -(1 - R/r) + r^2 cos(φ)^2 (dλ/dt)^2 
	... using dλ/dt = omega = angular velocity
	(ds/dt)^2 = -(1 - R/r) + r^2 cos(φ)^2 omega^2 
	... isn't ds^2 = -1 for unit objs in -+++ signature?
	-1/dt^2 = -(1 - R/r) + r^2 cos(φ)^2 omega^2 
	dt = 1 / sqrt(1 - R/r - r^2 cos(φ)^2 omega^2)
	... and i never remember if it's 1/ or not ...
	*/
	auto restingGravTimeDilations = map<real, real>(lats, [&](real lat){
		real r = geodeticDist(lat, 0, 0);
		return 1 / sqrt(1 - schwarzschildRadius_in_m / r);
	});
	
	auto sqr = [](auto x) { return x * x; };

	std::vector<real> gravTimeDilations;
	for (int i = 0; i < n; ++i) {
		auto lat = lats[i];
		auto beta = lorentzBetas[i];
		real r = geodeticDist(lat, 0, 0);
		real cosPhi = cos(degToRad(lat));
		gravTimeDilations.push_back(1 / sqrt(1 
			- schwarzschildRadius_in_m / r 
			- 1 / (1 - schwarzschildRadius_in_m / r) * beta * beta
		));
	}
	//print(matrix{lorentzGammas}:T())	-- verify how close to 1 all these are.  gnuplot just says '1' for everything.
	
	// alright, C++ here, so we can do the gnuplot in a second pass
	std::ofstream f("data.txt");
	f.precision(30);
	for (int i = 0; i < n; ++i) {
		f << i
			<< "\t" << lats[i]
			<< "\t" << dists[i]
			<< "\t" << vels[i]
			<< "\t" << lorentzGammas[i]
			<< "\t" << restingGravTimeDilations[i]
			<< "\t" << gravTimeDilations[i]
			<< std::endl;
	}

#if 0
	gnuplot{
		terminal = 'svg size 1024,768',
		output = 'dist.svg',
		style = 'data lines',
		data = {lats, dists},
		xrange = {-90, 90},
		xlabel = 'latitude',
		ylabel = 'meters',
		title='Earth surface distance from center, in meters',
		{using='1:2', title=''},
	}

	gnuplot{
		terminal = 'svg size 1024,768',
		output = 'vel.svg',
		style = 'data lines',
		data = {lats, vels},
		xrange = {-90, 90},
		xlabel = 'latitude',
		ylabel = 'meters per second',
		title='Earth surface speed about axis, in meters per second',
		{using='1:2', title=''},
	}

	gnuplot{
		terminal = 'svg size 1024,768',
		output = 'lorentzgamma.svg',
		style = 'data lines',
		data = {lats, lorentzGammas},
		xrange = {-90, 90},
		xlabel = 'latitude',
		ylabel = 'Lorentz gamma',
		format = {
			y = '%.20f',	-- can I label it 1 - %f?
		},
		title = 'Earth surface Lorentz boost',
		{using='1:2', title=''},
	}

	gnuplot{
		terminal = 'svg size 1024,768',
		output = 'grav-rest-dt.svg',
		style = 'data lines',
		data = {lats, restingGravTimeDilations},
		xrange = {-90, 90},
		xlabel = 'latitude',
		ylabel = 'dt',
		format = {
			y = '%.20f',	-- can I label it 1 - %f?
		},
		title = 'Earth surface gravitational time dilation (not considering Earth rotation)',
		{using='1:2', title=''},
	}

	gnuplot{
		terminal = 'svg size 1024,768',
		output = 'gravdt.svg',
		style = 'data lines',
		data = {lats, gravTimeDilations},
		xrange = {-90, 90},
		xlabel = 'latitude',
		ylabel = 'dt',
		format = {
			y = '%.20f',	-- can I label it 1 - %f?
		},
		title = 'Earth surface gravitational time dilation, considering rotation',
		{using='1:2', title=''},
	}
#endif
}
