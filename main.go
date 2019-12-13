package main

import "math"

// Variables are set of concetration of key components. Intended to be used as slice or array
type Variables struct {
	x float64 //x represents concentration of protonated PMO SM
	y float64 //y represents concentration of protonated base
	z float64 //z represents concentration of PMO product
	h float64 //h represents concentration of monomer hydrolyzed compound
}

//X is a slice of set of variables (concentrations). To be assigned in main()
var X []Variables

//Initial are set of initial concentration of key components
type Initial struct {
	PH0  float64 //PH0 represents initial concentration of PMO SM: M
	MCl0 float64 // MCl0 represents initial concentration of the monomer: M
	B0   float64 //B0 represents initial concentration of the base: M
}

//Ini is set of initial concentration of PMO, monomer, and base.
var Ini = Initial{
	PH0:  0.04363, // M
	MCl0: 0.06108, // M
	B0:   0.07329, //M
}

// Constants is set of kinetic and equilibrium constants
type Constants struct {
	KC float64 //KC is coupling reaction constant
	KS float64 //KS is monomer hydrolysis reaction
	K1 float64 //K1 is PMO SM protonation equilibrium constant
	K2 float64 //K2 is base protonation equilibrium constant
}

//Cons is a set of constants
var Cons = Constants{
	KC: 1.0,
	KS: 0.1,
	K1: 4.365e-9,
	K2: 1.7782e-11,
}

func main() {

}

func dzdt(Cons Constants, Ini Initial, X Variables) (difZ float64) {

	difZ = Cons.KC * (Ini.PH0 - X.x - X.z) * (Ini.MCl0 - X.z - X.h)

	return

}

func dhdt(Cons Constants, Ini Initial, X Variables) (difH float64) {
	difH = Cons.KS * (Ini.MCl0 - X.z - X.h)
	return
}

func equilibrium(Cons Constants, Ini Initial, X *Variables) {
	K := Cons.K1 / Cons.K2
	H := X.z + X.h
	R := Ini.PH0 - X.z
	B := Ini.B0 - X.z - X.h

	X.x = (-K*B - H - R + math.Sqrt(math.Pow(K*B+H+R, 2.0)+4.0*(K-1.0)*H*R)) /
		(2.0 * (K - 1.0))

	X.y = X.z + X.h - X.x
}
