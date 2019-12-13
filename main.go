package main

import (
	"fmt"
	"math"
)

// Variables are set of concetration of key components. Intended to be used as slice or array
type Variables struct {
	x float64 //x represents concentration of protonated PMO SM
	y float64 //y represents concentration of protonated base
	z float64 //z represents concentration of PMO product
	h float64 //h represents concentration of monomer hydrolyzed compound
}

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
	var div /*k,*/, iPC /*, back*/ int
	var tTot, dt float64

	tTot = 20.0
	div = 10000
	dt = tTot / float64(div)

	X := make([]Variables, div+1)

	iPC = 0
	X[0] = Variables{
		x: 0.0,
		y: 0.0,
		z: 0.0,
		h: 0.0,
	}

	for iPC = 1; iPC <= 3; iPC++ {
		var dif1, dif2, dif3, dif4 Variables
		var temp2, temp3, temp4 Variables

		dif1.z = dzdt(Cons, Ini, X[iPC-1])
		dif1.h = dhdt(Cons, Ini, X[iPC-1])
		temp2.z = X[iPC-1].z + dt*dif1.z/2.0
		temp2.h = X[iPC-1].h + dt*dif1.h/2.0
		equilibrium(Cons, Ini, &temp2)

		dif2.z = dzdt(Cons, Ini, temp2)
		dif2.h = dhdt(Cons, Ini, temp2)
		temp3.z = X[iPC-1].z + dt*dif2.z/2.0
		temp3.h = X[iPC-1].h + dt*dif2.h/2.0
		equilibrium(Cons, Ini, &temp3)

		dif3.z = dzdt(Cons, Ini, temp3)
		dif3.h = dhdt(Cons, Ini, temp3)
		temp4.z = X[iPC-1].z + dt*dif3.z
		temp4.h = X[iPC-1].h + dt*dif3.h
		equilibrium(Cons, Ini, &temp4)

		dif4.z = dzdt(Cons, Ini, temp4)
		dif4.h = dhdt(Cons, Ini, temp4)

		X[iPC].z = X[iPC-1].z + (dif1.z+dif2.z*2.0+dif3.z*2.0+dif4.z)*dt/6.0
		X[iPC].h = X[iPC-1].h + (dif1.h+dif2.h*2.0+dif3.h*2.0+dif4.h)*dt/6.0
		equilibrium(Cons, Ini, &X[iPC])
	}
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

	fmt.Println("ratio y/x:", X.y/X.x)
}
