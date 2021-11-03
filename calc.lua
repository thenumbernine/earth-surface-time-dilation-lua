#!/usr/bin/env luajit

local matrix = require 'matrix'
local gnuplot = require 'gnuplot'

local speedOfLight_in_m_per_s = 299792458
print("speed of light, in m/s:", speedOfLight_in_m_per_s)

local gravitationalConstant_in_m3_per_kg_s2 = 6.67384e-11
print("gravitational constant, in m^3 / (kg s^2):", gravitationalConstant_in_m3_per_kg_s2)


local sun = {
	name = 'Sun',
	equatorialRadius_in_m = 6.960e+8,
	inverseFlattening = 20000.0,
	--[[ from wikipedia ...
	https://en.wikipedia.org/wiki/Solar_rotation
	omega = 14.713 - 2.396 * sin(phi)^2 - 1.787 * sin(phi)^4 -- angular velocity in degrees per day, phi = latitude
	--]]
	rotationPeriod_in_s = 27 * 24 * 60 * 60,	-- varies by latitude really ...
	mass_in_kg = 1.9891e+30,
}
local earth = {
	name = 'Earth',
	equatorialRadius_in_m = 6378136,
	inverseFlattening = 298.257223563,
	rotationPeriod_in_s = (23 + (56 + 4.09053083288 / 60) / 60) * 60 * 60,	-- sidereeal?
	mass_in_kg = 5.9736e+24,
}
local mercury = {
	name = 'Mercury',
	equatorialRadius_in_m = 2440e+3,
	inverseFlattening = math.huge,
	rotationPeriod_in_s = 58.6462 * 60 * 60 * 24,
	mass_in_kg = 3.302e+23,
}
local venus = {
	name = 'Venus',
	equatorialRadius_in_m = 6051.893e+3,
	inverseFlattening = math.huge,
	rotationPeriod_in_s = -243.0185 * 60 * 60 * 24,
	mass_in_kg = 4.8685e+24,
}
local moon = {
	name = 'Moon',
	equatorialRadius_in_m = 1.73753e+6,
	inverseFlattening = 1/0.0012,
	mass_in_kg = 7.349e+22,
	rotationPeriod_in_s = 27.321662 * 60 * 60 * 24,
}
local mars = {
	name = 'Mars',
	equatorialRadius_in_m = 3397e+3,
	inverseFlattening = 154.409,
	mass_in_kg = 6.4185e+23,
	rotationPeriod_in_s = 24.622962/24 * 60 * 60 * 24,
}
local jupiter = {
	name = 'Jupiter',
	equatorialRadius_in_m = 71492e+3,
	inverseFlattening = 1/0.06487,
	mass_in_kg = 1.89813e+27,
	rotationPeriod_in_s = 10 * 60 * 60,	-- 10 hours
}
local saturn = {
	name = 'Saturn',
	equatorialRadius_in_m = 60268e+3,
	inverseFlattening = 1/0.09796,
	mass_in_kg = 5.68319e+26,
	rotationPeriod_in_s = 24 + 60 * (39 + 60 * 10),	-- https://www.jpl.nasa.gov/news/rotation-period-of-saturn-determined
}
local uranus = {
	name = 'Uranus',
	equatorialRadius_in_m = 25559e+3,
	inverseFlattening = 1/0.02293,
	mass_in_kg = 8.68103e+25,
	rotationPeriod_in_s = 23 * 60 * 60,	-- 23 hours (+5 / -2 hours ... ?)
}
local neptune = {
	name = 'Neptune',
	equatorialRadius_in_m = 24766e+3,
	inverseFlattening = 1/0.0171,
	mass_in_kg = 1.0241e+26,
}
local pluto = {
	name = 'Pluto',
	equatorialRadius_in_m = 1.151e+6,
	rotationPeriod_in_s = 6.387230 * 60 * 60 * 24,	-- but it's negative
	mass_in_kg = 1.314e+22,
	inverseFlattening = math.huge,
}
local FakePlanet = {
	name = 'Fake',
	equatorialRadius_in_m = 1e+4,
	inverseFlattening = 1e+2,
	rotationPeriod_in_s = 60,
	mass_in_kg = 1e+30,
}

--local planet = sun
local planet = earth
--local planet = saturn



print('planet', planet.name)

-- check for earth
local equatorialRadius_in_m = planet.equatorialRadius_in_m
print("equitorialRadius_in_m:", planet.equatorialRadius_in_m)

-- check for earth for sidereal period
local rotationPeriod_in_s = planet.rotationPeriod_in_s
print("rotationPeriod_in_s:", planet.rotationPeriod_in_s)

-- check for earth
local mass_in_kg = planet.mass_in_kg
print('mass_in_kg:', planet.mass_in_kg)

-- check for earth
local schwarzschildRadius_in_m = 2 * mass_in_kg * gravitationalConstant_in_m3_per_kg_s2 / (speedOfLight_in_m_per_s * speedOfLight_in_m_per_s)
print('schwarzschildRadius_in_m:', schwarzschildRadius_in_m)

-- check for earth
local inverseFlattening = planet.inverseFlattening
print("inverse flattening:", inverseFlattening)

local eccentricitySquared = (2 * inverseFlattening - 1) / (inverseFlattening * inverseFlattening)
print('eccentricitySquared', eccentricitySquared)

local angularSpeed_zAxis_in_1_per_s = (2 * math.pi) / rotationPeriod_in_s
print('angularSpeed_zAxis_in_1_per_s', angularSpeed_zAxis_in_1_per_s)

local rotationAxis = matrix{0,0,1}
print('rotationAxis', rotationAxis)

local angularVelocity_in_1_per_s = rotationAxis * angularSpeed_zAxis_in_1_per_s
print('angularVelocity_in_1_per_s', angularVelocity_in_1_per_s)

-- ellipsoid inertia matrix
local a = equatorialRadius_in_m
local b = a
local c = a * (eccentricitySquared * inverseFlattening - 1)
local Inertia_in_kg_m2 = matrix{
	{mass_in_kg * (b*b + c*c) / 5, 0, 0},
	{0, mass_in_kg * (a*a + c*c) / 5, 0},
	{0, 0, mass_in_kg * (a*a + b*b) / 5},
}
-- check for earth
print('Inertia_in_kg_m2')
print(Inertia_in_kg_m2)

-- check for earth
local angularMomentum_in_kg_m2_per_s = Inertia_in_kg_m2 * angularVelocity_in_1_per_s
print('angularMomentum_in_kg_m2_per_s', angularMomentum_in_kg_m2_per_s)

-- ~4, for the earth.  who would've thought.
local kerr_a_in_m = rotationAxis * angularMomentum_in_kg_m2_per_s / mass_in_kg / speedOfLight_in_m_per_s
print('kerr_a_in_m', kerr_a_in_m)
-- J / (M c) = a << r is a requirement for the kerr metric to work

-- lat lon in degrees
local function geodeticPos(lat, lon, height)
	local phi = math.rad(lat)
	local lambda = math.rad(lon)
	local cosPhi = math.cos(phi)
	local sinPhi = math.sin(phi)
	local sinPhiSquared = sinPhi * sinPhi
	local N = equatorialRadius_in_m / math.sqrt(1 - eccentricitySquared * sinPhiSquared)
	local NPlusH = N + height
	local x = NPlusH * cosPhi * math.cos(lambda)
	local y = NPlusH * cosPhi * math.sin(lambda)
	local z = (N * (1 - eccentricitySquared) + height) * sinPhi
	return x, y, z	-- in meters
end

--[[
is this the isotropic r or the Schwarzschild r? 

from "catalog of spacetimes" https://arxiv.org/pdf/0904.4184.pdf
sch_r = iso_r (1 + r_s / (4 iso_r))^2

therefore:
iso_r^2 + iso_r (r_s/2 - sch_r) + r_s^2/16 = 0
iso_r = (sch_r - r_s/2 +- sqrt(sch_r (sch_r - r_s)))/2
so that when r_s ~ 0, iso_r = sch_r

alright, catalog of spacetimes, 
2.2.2 pseudo-Cartesian: sch_r^2 = x^2 + y^2 + z^2
2.2.32 isotropic coordiantes: iso_r^2 = x^2 + y^2 + z^2
so ... we are dealing with different sets of (x,y,z)'s for different metrics.
so which one is our Newtonian (x,y,z)?

from wiki on istropic coordinates "This means that ... nor does the radial coordinate faithfully represent radial distances"

Because isotropic does have its own spherical and cartesian form, 
and Schwarzschild non-isotropic also has its own spherical and (pseudo)cartesian form,
I'm guessing that the Newtonian 'r' is the same as the Schwarzschild 'r', and not the isotropic 'r' ... ???
--]]
local function geodeticDist(lat, lon, height)
	local x,y,z = geodeticPos(lat, lon, height)
	return math.sqrt(x*x + y*y + z*z)
end

local function geodeticVel(lat, lon, height)
	local phi = math.rad(lat)
	local lambda = math.rad(lon)
	local cosPhi = math.cos(phi)
	local sinPhi = math.sin(phi)
	local sinPhiSquared = sinPhi * sinPhi
	local N = equatorialRadius_in_m / math.sqrt(1 - eccentricitySquared * sinPhiSquared)
	local NPlusH = N + height
	local x = -angularSpeed_zAxis_in_1_per_s * NPlusH * cosPhi * math.sin(lambda)
	local y = angularSpeed_zAxis_in_1_per_s * NPlusH * cosPhi * math.cos(lambda)
	local z = 0
	return x, y, z
end

local function geodeticSpeed(lat, lon, height)
	local x, y, z = geodeticVel(lat, lon, height)
	return math.sqrt(x*x + y*y + z*z)
end

local n = 180
local x = matrix{n}:lambda(function(i) return (i-1)/(n-1) end)	-- on edge
local lats = x * 180 - 90

local function plotvec(args)
	local col = args[1]
	assert(matrix:isa(col))
	args.terminal = 'svg size 1024,768 background rgb "white"'
	args.style = 'data lines'
	args.data = {lats, col}
	args.xlabel = 'latitude (degrees)'
	args.xrange = {-90, 90}
	args.title = planet.name..' '..args.title
	args[1] = {using='1:2', title=''}
	gnuplot(args)
end

local dists = lats:map(function(lat) return geodeticDist(lat, 0, 0) end)
-- check
plotvec{dists, output='dist.svg', ylabel='distance (meters)', title='surface distance from center'}

local vels = lats:map(function(lat) return geodeticSpeed(lat, 0, 0) end)
-- check
plotvec{vels, output='vel.svg', ylabel='speed (meters per second)', title='surface speed about axis'}

-- special relativity beta and gamma ...
local lorentzBetas = vels / speedOfLight_in_m_per_s
-- check
plotvec{lorentzBetas, output='lorentzbeta.svg', ylabel='Lorentz beta', format={y='%.20f'}, title='surface Lorentz beta'}

local lorentzGammas = lorentzBetas:map(function(beta) return math.sqrt(1 - beta * beta) end)
-- check
plotvec{lorentzGammas, output='lorentzgamma.svg', ylabel='Lorentz gamma', format={y='%.20f'}, title='surface Lorentz gamma'}

local specialRelativisticChangeInDayInNS = (lorentzGammas - 1) * 24 * 60 * 60 * 1e+9
-- check
plotvec{specialRelativisticChangeInDayInNS, output='change-in-day-sr.svg', ylabel='dt (ns)', title='change in duration of a day'}


--[[
GR 4-vel ... how to calculate ... seems I remember seeing a few defs in the GRHD papers ...
has to be such that the relativistic norm is -1 ...
ds^2 = (-alpha^2 + beta_i beta^i) u^0 + 2 gamma_ij beta^i u^j + gamma_ij u^i u^j
-1 = g_00 (u^0)^2 + 2 g_0i u^i u^0 + g_ij u^i u^j

options ...
1) pick u^i, calculate u^0
g_00 (u^0)^2 + 2 g_0i u^i u^0 + g_ij u^i u^j + 1 = 0
u^0 = 1/(2 * g_00) (
	-2 g_0i u^i
	+- sqrt(
		4 (g_0i u^i)^2
		- 4 * g_00 * (g_ij u^i u^j + 1)
	)
)

2) pick u^0, calculate |u^i|
g_00 (u^0)^2 + 1 + 2 g_0i uhat^i u^0 |u^i| + g_ij uhat^i uhat^j |u^i|^2 = 0
|u^i| = 1/(2 * (g_00 (u^0)^2 + 1)) (
	-2 g_0i uhat^i u^0
	+- sqrt(
		4 (g_0i uhat^i u^0)^2
		- 4 (g_00 (u^0)^2 + 1) (g_ij uhat^i uhat^j)
	)
)

I like option 1
SR's special-case of #1 would be considering u^i to be v^i/c, can I use the same trick in GR?
but the downside to option #1 is ... should the spatial hypersurface metric influence the 3-norm of the velocity?


... so which is it?

for an object at rest, this is easy:
u^0 = sqrt(-1 / g_00)

how about an object in motion?
I'd say to use u^i = beta^i ... but for luminal velocities, i.e. null geodesics ... 
0 = g_00 (u^0)^2 + 2 g_0i u^0 u^j + g_ij u^i u^j
... in a diagonal metric ...
(u^0)^2 = g_ii (u^i)^2 / -g_00
u^0 = sqrt( sum_i g_ii (u^i)^2 / -g_00 )
... but does this mean that u^i can reach the value of 1?

what about for a timelike geodesic in a diagonal metric?
-1 = g_00 (u^0)^2 + g_ii (u^i)^2
--]]

local function newtonApproxMetric(x,y,z)
	local r = math.sqrt(x*x + y*y + z*z)
	-- would this 'r' be isotropic or not?
	local R = schwarzschildRadius_in_m
	local neg2phi = R / r 
	local g00 = -1 + neg2phi
	local gii = 1 + neg2phi
	return matrix{
		{g00, 0, 0, 0},
		{0, gii, 0, 0},
		{0, 0, gii, 0},
		{0, 0, 0, gii},
	}
end

local function schwarzschildSphericalMetric(x,y,z)
	local theta = math.atan2(math.sqrt(x*x + y*y), z)
	-- is this isotropic r or schwarzschild r?
	local r = math.sqrt(x*x + y*y + z*z)
	--[[ if it's isotropic r then ...
	r = r * (1 + R / (4 * r))^2
	--]]
	local R = schwarzschildRadius_in_m
	local rsinth = r * math.sin(theta)
	local ir = 1 / r
	local u = 1 - R * ir
	return matrix{
		{-u, 0, 0, 0},
		{0, 1/u, 0, 0},
		{0, 0, r*r, 0},
		{0, 0, 0, rsinth*rsinth},
	}
end

-- if non-isotropic is pseudo-cartesian, then does that mean isotropic is genuine cartesian?
local function schwarzschildPseudoCartesianMetric(x,y,z)
	local r = math.sqrt(x*x + y*y + z*z)
	--[[ if |xyz|2 is isotropic then ...
	r = r * (1 + R / (4 * r))^2
	--]]
	local R = schwarzschildRadius_in_m
	local ir = 1 / r
	local ux = x * ir
	local uy = y * ir
	local uz = z * ir
	local irsq = ir * ir
	local u = 1 - R * ir		-- what exactly is this value?  sqrt(1 - (escape vel / speed of light)^2) = ... what?
	local iu = 1 / u
	local Rsq = R * R
	local g00 = -u
	local uxsq = ux * ux
	local uysq = uy * uy
	local uzsq = uz * uz
	local g11 = uxsq * iu + uysq + uzsq
	local g22 = uxsq + uysq * iu + uzsq
	local g33 = uxsq + uysq + uzsq * iu
	local g12 = ux * uy * iu * Rsq * ir
	local g23 = uy * uz * iu * Rsq * ir
	local g31 = uz * ux * iu * Rsq * ir
	return matrix{
		{g00, 0, 0, 0},
		{0, g11, g12, g31},
		{0, g12, g22, g23},
		{0, g31, g23, g33},
	}
end

local function isotropicCartesianMetric(x,y,z)
	local r = math.sqrt(x*x + y*y + z*z)
	local R = schwarzschildRadius_in_m
	-- [[ if it's schwarzschild r, then convert it to isotropic r:
	r = .5 * (r - .5 * R + math.sqrt(r * (r - R)))
	--]]
	local ir = 1 / r
	local Rr = .25 * R * ir
	local mup = 1 + Rr
	local mum = 1 - Rr
	local mup2 = mup * mup 
	local mum2 = mum * mum
	local mup4 = mup2 * mup2
	return matrix{
		{-mum2/mup2, 0, 0, 0},
		{0, mup4, 0, 0},
		{0, 0, mup4, 0},
		{0, 0, 0, mup4},
	}
end

local function isotropicSphericalMetric(x,y,z)
	local r = math.sqrt(x*x + y*y + z*z)
	local theta = math.atan2(math.sqrt(x*x + y*y), z)	-- is theta pre or post isotropic radial transform?
	local R = schwarzschildRadius_in_m
	-- [[ if it's schwarzschild r, then convert it to isotropic r:
	r = .5 * (r - .5 * R + math.sqrt(r * (r - R)))
	--]]
	local ir = 1 / r
	local Rr = .25 * R * ir
	local mup = 1 + Rr
	local mum = 1 - Rr
	local mup2 = mup * mup 
	local mum2 = mum * mum
	local mup4 = mup2 * mup2
	local rsinth = r * math.sin(theta)	-- spherical sin(theta) = geo cos(lat)
	return matrix{
		{-mum2/mup2, 0, 0, 0},
		{0, mup4, 0, 0},
		{0, 0, mup4*r*r, 0},
		{0, 0, 0, mup4*rsinth*rsinth},
	}
end

local function kerrBoyerLindquistSphericalMetric(x,y,z)
	local r = math.sqrt(x*x + y*y + z*z)
	local theta = math.atan2(math.sqrt(x*x + y*y), z)
	local R = schwarzschildRadius_in_m
	local rsq = r * r
	local a = kerr_a_in_m
	local asq = a * a
	local sinth = math.sin(theta)
	local costh = math.cos(theta)
	local sinthsq = sinth * sinth
	local costhsq = costh * costh
	local sigma = rsq + asq * costhsq 
	local delta = rsq - R * r + asq
	local g00 = -(1 - R * r / sigma)
	local g02 = -R * a * r * sinthsq / sigma
	local g11 = sigma / delta
	local g22 = sigma
	local g33 = (rsq + asq + R * asq * r * sinthsq / sigma) * sinthsq
	return matrix{
		{g00, 0, g02, 0},
		{0, g11, 0, 0},
		{g02, 0, g22, 0},
		{0, 0, 0, g33},
	}
end

local eta = matrix{4,4}:lambda(function(i,j)
	return i == j and (i == 1 and -1 or 1) or 0
end)
local function kerrSchildMetric(x,y,z)
	local R = schwarzschildRadius_in_m
	local a = kerr_a_in_m
	local asq = a * a
	--[[
	r is defined by (x^2 + y^2) / (r^2 + a^2) + z^2 / r^2 = 1
	r^4 + r^2 (a^2 - (x^2 + y^2 + z^2)) - z^2 a^2 = 0
	r^2 = 1/2 (
		- (a^2 - (x^2 + y^2 + z^2))
		+- sqrt(
			(a^2 - (x^2 + y^2 + z^2))^2
			- 4 (- z^2 a^2)
		)
	)
	r^2 = 1/2 (
		(x^2 + y^2 + z^2) - a^2
		+- sqrt(
			(a^2 - (x^2 + y^2 + z^2))^2
			+ 4 z^2 a^2
		)
	)
	--]]
	local zsq = z*z
	local b = asq - (x*x + y*y + zsq)
	local rsq = .5 * (
		-b
		+ math.sqrt(
			b * b
			+ 4 * asq * zsq
		)
	)
	local r = math.sqrt(rsq)
	local rcub = r * rsq
	local rqu = rsq * rsq

	-- f = (2 G M r^3) / (r^4 + a^2 z^2)
	local f = R * rcub / (rqu + asq * zsq)

	-- k = (1, (r x + a y) / (r^2 + a^2), (r y - a x) / (r^2 + a^2), z / r)
	local invdenom = 1 / (rsq + asq)
	local k = matrix{1, (r * x + a * y) * invdenom, (r * y - a * x) * invdenom, z / r}

	-- g_uv = eta_uv + f k_u k_v
	local g = eta + f * matrix.outer(k, k)
	return g
end


--[[
n = u * G * u
to make n = -1 ...

u => u / sqrt(-uGu)
such that 
	u / sqrt(-uGu) * G * u / sqrt(-uGu)
	= u * G * u / (-uGu)
	= -1

however why isn't this always the best option?
because we don't always know what u^0 should be.
this is fine for purely-timelike vectors
but for spatial vectors, how to determine what the timelike component is?
--]]
--[[
local function normalize(u, g)
	local n = u * g * u
	return u / math.sqrt(math.abs(n))
end
--]]

--[[
normalize first, and then solve

-1 = g_00 (u^0)^2 + 2 g_0i u^0 u^i + g_ij u^i u^j

g_00 (u^0)^2 + 2 g_0i u^i u^0 + g_ij u^i u^j + 1 = 0

u^0 = (
	-g_0i u^i
	+- sqrt(
		(g_0i u^i)^2
		- g_00 * (g_ij u^i u^j + 1)
	)
) / g_00

for u^i = 0 ...

u^0 = +- sqrt(-g_00) / g_00
u^0 = -1/sqrt(-g_00)
looks like minus is our choice 
--]]
-- [[
local function normalize(u, g)
	u = matrix(u)
	u[1] = 0
	local g00 = g[1][1]
	local g0i_ui = g[1] * u
	local u3norm = u * g * u
	u[1] = (-g0i_ui - math.sqrt(g0i_ui * g0i_ui - g00 * (u3norm + 1))) / g00
	return u
end
--]]

--[[ for resting objects ...
Schwarzschild metric ...
	ds^2 = -(1 - R/r) dt^2 + 1/(1 - R/r) dr^2 + r^2 (dθ^2 + sin(θ)^2 dφ^2)
... at rest ...
	ds^2 = -(1 - R/r) dt^2
... assume ds^2 = -1 ...
	-1 = -(1 - R/r) dt^2
... isn't it bad math to assume ds^2 = -1?
wouldn't a better explanation be solving for -1 = |dx/ds|?
	-1 = -(1 - R/r) dt/ds^2
	dt/ds = 1 / sqrt(1 - R/r)
so this is our u^0 = dt/ds ... but is this our clock duration?
is the clock dt/ds = u^0 or ds/dt = 1/u^0?
	ds/dt = sqrt(1 - R/r)
u^0 = how fast the moving object travels through the time coordinate
whereas, when measuring a unit of resting time coordinate, 1/u^0 will be how much the moving clock ticks in a single unit of resting time coordinate ... 
... right?
that's why we use 1/u^0?
--]]
local metrics = {
	{
		g = newtonApproxMetric,
		title = 'Newtonian approximation',
	},
	{
		g = schwarzschildSphericalMetric,
		title = 'Schwarzschild spherical',
	},
	{
		g = schwarzschildPseudoCartesianMetric,
		title = 'Schwarzschild pseudo-cartesian',
	},
	
	-- how does ds/dt for resting objects look under isotropic coordinates?
	{
		g = isotropicCartesianMetric,
		title = 'isotropic Cartesian',
	},
	{
		g = isotropicSphericalMetric,
		title = 'isotropic spherical',
	},
	
	
	-- though, technically, a stationary point in Kerr metric would be equivalent to a rotating point on a Schwarzschild metric, right?
	{
		g = kerrBoyerLindquistSphericalMetric,
		title = 'Kerr Boyer-Lindquist spherical',
	},

	{
		g = kerrSchildMetric,
		title = 'Kerr-Schild',
	},
} 
-- TODO these are all nearly identical. just use one.
-- maybe aftewards plot the matrix of differences, but they are all nearly zero anyways.
for _,metric in ipairs(metrics) do
	local g = metric.g
	local title = metric.title
	local dsdts = lats:map(function(lat)
		local x,y,z = geodeticPos(lat, 0, 0)
		-- TODO instead of {1,0,0,0}, each metric should have its own associated timelike vector, equal to d/dt of the chart
		return 1 /  normalize(matrix{1,0,0,0}, g(x,y,z))[1]
	end)
	metric.dsdts = dsdts
	plotvec{
		dsdts,
		output = 'gr-rest-obj-dsdt-'..title..'.svg',
		ylabel = 'dt',
		format = {y = '%.20f'},
		title = 'surface gravitational time dilation of an object at rest (4-vel=[1,0,0,0]), '..title..' metric',
	}

--	local dsdtForDist = (x * 2 * equatorialRadius_in_m):map(function(r)
--	end)

	-- would the GR standard time be the pole or the equator?
	local dt00ref = dsdts[math.ceil(#dsdts/2)]
	--[[
	in order to normalize our day length 
	(since the schwarzschild metric g_00 term will always be non-unit on the surface at all latitudes
	i will choose the dt at the pole (which is stationary) as the reference for how long a day should last
	
	honestly if the previous graph is identical for each metric then this will be identical too
	--]]
	local changeInDayInNS = (dsdts / dt00ref - 1) * 24 * 60 * 60 * 1e+9
	metric.changeInDayInNS = changeInDayInNS
	plotvec{
		changeInDayInNS,
		output = 'change-in-day-gr-obj-at-rest-'..title..'.svg',
		ylabel = 'dt (ns)',
		title = 'change in day duration of object at rest, '..title..' metric',
	}

	
end
print('differences in metrics:')
print(matrix{#metrics, #metrics}:lambda(function(i,j)
	local mi = metrics[i]
	local mj = metrics[j]
	return (mi.dsdts - mj.dsdts):normSq()
end))


--[[ for rotating objects, in spherical coordinates
Schwarzschild metric ...
	ds^2 = -(1 - R/r) dt^2 + 1/(1 - R/r) dr^2 + r^2 (dθ^2 + sin(θ)^2 dφ^2)
... for an object not moving radially, not moving latitudally ...
	ds^2 = -(1 - R/r) dt^2 + r^2 sin(θ)^2 dφ^2 
... normalized to -1 in -+++ signature
	-1 = -(1 - R/r) (dt/ds)^2 + r^2 sin(θ)^2 (dφ/ds)^2 
	-(ds/dt)^2 = -(1 - R/r) + r^2 sin(θ)^2 (dφ/dt)^2 
	-(ds/dt)^2 = -(1 - R/r) + (r sin(θ) dφ/dt)^2 
.. substitute dφ/dt = ω, solve for dt
	-(ds/dt)^2 = -(1 - R/r) + (r sin(θ) ω)^2 
	(ds/dt)^2 = 1 - R/r - (r sin(θ) ω)^2
	ds/dt = sqrt(1 - R/r - (r sin(θ) ω)^2)
... but now we have some weird speed of light constraint on our ω ...

	is beta = v / c = r sin(θ) ω / c ?
	not necessarily, since at equator, r = (N + h), at pole r = (N (1 - e^2) + h) ... so ...
	... should I be doing this in pseudo-Cartesian coordinates?
	or like the paper here: https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=904814
	... doing it in isotropic coordinates?

I guess I could try to claim that r dφ/dt needs to be converted to a 4-velocity ...
I'm not convinced on the schwarzschild implementation, esp with however the 4-velocity should be converted ...
--]]
local dt00ref = math.sqrt(1 - schwarzschildRadius_in_m / geodeticDist(90, 0, 0))
local gravTimeDilationSchwarzschildSpherical = matrix{n}:lambda(function(i)
	local lat = lats[i]
	local beta = lorentzBetas[i]
	local x,y,z = geodeticPos(lat, 0, 0)
	--[[
	local cosPhi = math.cos(math.rad(lat))
	local r = math.sqrt(x*x + y*y + z*z)
	--]]
	local g = schwarzschildSphericalMetric(x,y,z)
	return math.sqrt(
		-g[1][1]
		-- how did I come up with part this again?
		- beta * beta
	) / dt00ref
end)
gnuplot{
	terminal = 'svg size 1024,768',
	output = 'grav-dt-schwarzschild-spherical.svg',
	style = 'data lines',
	data = {lats, gravTimeDilationSchwarzschildSpherical},
	xrange = {-90, 90},
	xlabel = 'latitude',
	ylabel = 'dt',
	format = {
		y = '%.20f',	-- can I label it 1 - %f?
	},
	title = planet.name..' surface gravitational time dilation, considering rotation, Schwarzschild',
	{using='1:2', title=''},
}

--[[
ds^2 = g_uv dx^u dx^v
-1 = g_uv dx^u/ds dx^v/ds
-1 = g_00 (dt/ds)^2 + 2 g_0i dt/ds dx^i/ds + g_ij dx^i/ds dx^j/ds

for schwarzschild pseudo-carteisan, g_0i = 0
dt/ds = u^0 = sqrt( (g_ij dx^i/ds dx^j/ds + 1) / -g_00 )
ds/dt = sqrt( -g_00 / (g_ij dx^i/ds dx^j/ds + 1) )
--]]
local gravTimeDilationSchwarzschildPseudoCartesian = matrix{n}:lambda(function(i)
	local lat = lats[i]
	local x,y,z = geodeticPos(lat,0,0)
	local vx,vy,vz = geodeticVel(lat,0,0)
	local c = speedOfLight_in_m_per_s
	local u3 = matrix{0, vx/c, vy/c, vz/c}			-- is this right ...
	local g = schwarzschildPseudoCartesianMetric(x,y,z)
	local ds3 = u3 * g * u3
	local u0 = math.sqrt((ds3 + 1) / -g[1][1])
	--[[ such that now u0,u3 should give us -1 everywhere ... check
	local u = matrix{u0, u3[2], u3[3], u3[4]}
	return u * g * u
	--]]
	-- so why is this so out of whack?
--	error'here'
	return (1/u0) / dt00ref
end)
plotvec{
	gravTimeDilationSchwarzschildPseudoCartesian,
	output = 'grav-dt-schwarzschild-pseudocartesian.svg',
	ylabel = 'dt',
	format = {y='%.20f'},
	title = 'surface gravitational time dilation, considering rotation, Schwarzschild pseudo-Cartesian',
}
