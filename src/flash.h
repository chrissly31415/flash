/*
 * flash.h
 *
 *  Created on: Jan 6, 2011
 *  Author: Christoph Loschen
 */

#ifndef FLASH_H_
#define FLASH_H_

//gas constant [kcal/mol]
#define R 0.0019858775
//z=0 is the interface -> z_max+1
#define zmax 11
//monomer length in nanometer
//4,4-MDI ~ 13 Angstrom = 13 * 10-10m = 1.3 nm
//PEG O-C-C-O 3.1A cis and 3.8 trans position
#define mlength 1.0

//hexagonal lattice
#define cubic1 0.166666667
#define cubic2 0.666666667

class flash {
public:
	flash();
	virtual ~flash();
	//static variables
};


#endif /* FLASH_H_ */
