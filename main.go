package main

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
	kc float64 //kc is coupling reaction constant
	ks float64 //ks is monomer hydrolysis reaction
	K1 float64 //K1 is PMO SM protonation equilibrium constant
	K2 float64 //K2 is base protonation equilibrium constant
}

func main() {

}