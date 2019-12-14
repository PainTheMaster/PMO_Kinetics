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
	KC: 38.93, // conversion at 5 hours is 99% when KC == 38.93
	KS: 0.0,
	K1: 4.365e-9,
	K2: 1.7782e-11,
}

func main() {
	var div /*k,*/, iPC /*, back*/ int
	var tTot, dt float64

	tTot = 20.0
	div = 1000
	dt = tTot / float64(div)

	X := make([]Variables, div+1)

	iPC = 0
	X[0] = Variables{
		x: 0.0,
		y: 0.0,
		z: 0.0,
		h: 0.0,
	}

	for iPC = 1; iPC <= div; iPC++ {
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

	fmt.Println("RK Done!")

	Y := make([]Variables, div+1)

	for iPC = 0; iPC <= 3; iPC++ {
		Y[iPC] = X[iPC]
	}

	for iPC = 4; iPC <= div; iPC++ {
		//Prediction
		predZ0 := dzdt(Cons, Ini, Y[iPC-4])
		predZ1 := dzdt(Cons, Ini, Y[iPC-3])
		predZ2 := dzdt(Cons, Ini, Y[iPC-2])
		predZ3 := dzdt(Cons, Ini, Y[iPC-1])

		predH0 := dhdt(Cons, Ini, Y[iPC-4])
		predH1 := dhdt(Cons, Ini, Y[iPC-3])
		predH2 := dhdt(Cons, Ini, Y[iPC-2])
		predH3 := dhdt(Cons, Ini, Y[iPC-1])

		Y[iPC].z = Y[iPC-1].z + (55.0*predZ3-59.0*predZ2+37.0*predZ1-9.0*predZ0)*dt/24.0
		Y[iPC].h = Y[iPC-1].h + (55.0*predH3-59.0*predH2+37.0*predH1-9.0*predH0)*dt/24.0

		equilibrium(Cons, Ini, &Y[iPC])

		//Correction
		corrZ1 := dzdt(Cons, Ini, Y[iPC-3])
		corrZ2 := dzdt(Cons, Ini, Y[iPC-2])
		corrZ3 := dzdt(Cons, Ini, Y[iPC-1])
		corrZ4 := dzdt(Cons, Ini, Y[iPC-0])

		corrH1 := dhdt(Cons, Ini, Y[iPC-3])
		corrH2 := dhdt(Cons, Ini, Y[iPC-2])
		corrH3 := dhdt(Cons, Ini, Y[iPC-1])
		corrH4 := dhdt(Cons, Ini, Y[iPC-0])

		Y[iPC].z = Y[iPC-1].z + (9.0*corrZ4+19.0*corrZ3-5.0*corrZ2+corrZ1)*dt/24.0
		Y[iPC].h = Y[iPC-1].h + (9.0*corrH4+19.0*corrH3-5.0*corrH2+corrH1)*dt/24.0

		equilibrium(Cons, Ini, &Y[iPC])
	}

	for i := range X {
		fmt.Printf("%f,%.10f,%.10f\n", float64(i)*dt, Ini.PH0-X[i].z, Ini.PH0-Y[i].z)
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
}
