import sys
import json
from collections import OrderedDict

class Parameters():
    """Container for forcefield parameters
    """
    def __init__(self):

        self.atom_params = OrderedDict()
        self.atom_params['atoms'] = []
    
        self.bond_params = OrderedDict()
        self.bond_params['bonds'] = []
        
        self.angle_params = OrderedDict()
        self.angle_params['angles'] = []

        self.dihedral_params = OrderedDict()
        self.dihedral_params['dihedrals'] = []

        self.improper_params = OrderedDict()
        self.improper_params['impropers'] = []
        
        
        self.n_params = OrderedDict()
        self.n_params['n_params'] = []
        self.n_params['n_params'].append(OrderedDict([('n_atoms', 0),
                                                      ('n_bonds', 0),
                                                      ('n_angles', 0),
                                                      ('n_dihedrals', 0),
                                                      ('n_impropers', 0)]))

        self.n_atoms = 0
        self.n_bonds = 0
        self.n_angles = 0
        self.n_dihedrals = 0
        self.n_impropers = 0

    def add_atom(self, name, bond_type, atomic_number, mass, charge, ptype, sigma, epsilon):
        self.n_atoms += 1
        self.n_params['n_params'][0]['n_atoms']= self.n_atoms;
        
        self.atom_params['atoms'].append(OrderedDict([
             ('name', name),
             ('bond_type', bond_type),
             ('atomic_number', atomic_number),
             ('mass', mass),
             ('charge', charge),
             ('ptype', ptype),
             ('sigma', sigma),
             ('epsilon', epsilon)
        ]))
    
    def add_bond_harmonic(self, i, j, b0, kb):
        self.n_bonds += 1
        self.n_params['n_params'][0]['n_bonds']= self.n_bonds

        self.bond_params['bonds'].append(OrderedDict([
            ('type', 'harmonic'),
            ('i', i),
            ('j', j),
            ('func', 1),
            ('b0', b0),
            ('kb', kb)
        ]))
   
    def add_bond_G96(self, i, j, b0, kb):
        self.n_bonds += 1
        self.n_params['n_params'][0]['n_bonds']= self.n_bonds
        
        self.bond_params['bonds'].append(OrderedDict([
                                         ('type', 'G96'),
                                         ('i', i),
                                         ('j', j),
                                         ('func', 2),
                                         ('b0', b0),
                                         ('kb', kb)
                                         ]))
    def add_bond_fene(self, i, j, bm, kb):
        self.n_bonds += 1
        self.n_params['n_params'][0]['n_bonds']= self.n_bonds
        
        self.bond_params['bonds'].append(OrderedDict([
                                         ('type', 'fene'),
                                         ('i', i),
                                         ('j', j),
                                         ('func', 7),
                                         ('bm', bm),
                                         ('kb', kb)
                                         ]))
    def add_bond_morse(self, i, j, b0, D, beta):
        self.n_bonds += 1
        self.n_params['n_params'][0]['n_bonds']= self.n_bonds
        
        self.bond_params['bonds'].append(OrderedDict([
                                         ('type', 'morse'),
                                         ('i', i),
                                         ('j', j),
                                         ('func', 3),
                                         ('b0', b0),
                                         ('D', D),
                                         ('beta', beta)
                                         ]))

    def add_angle_harmonic(self, i, j, k, theta_0, kb):
        self.n_angles += 1
        self.n_params['n_params'][0]['n_angles']= self.n_angles;
        
        self.angle_params['angles'].append(OrderedDict([
            ('type', 'harmonic'),
            ('i', i),
            ('j', j),
            ('k', k),
            ('func', 1),
            ('theta_0', theta_0),
            ('kb', kb)
         ]))
    
    def add_angle_G96(self, i, j, k, theta_0, kb):
         self.n_angles += 1
         self.n_params['n_params'][0]['n_angles']= self.n_angles;
         
         self.angle_params['angles'].append(OrderedDict([
                                            ('type', 'G96'),
                                            ('i', i),
                                            ('j', j),
                                            ('k', k),
                                            ('func', 2),
                                            ('theta_0', theta_0),
                                            ('kb', kb)
                                            ]))
    '''Gromacs type 3 (i.e., func =3) Ryckaert-Bellemans (RB) dihedrals
        '''
    #to convert OPLS to RB:
    #OPLS form = 0.5*(F1*(1+cos(phi))+ F2*(1-cos(2phi))+ F3*(1+cos(3phi))+F4*(1-cos(4phi)))
    #RB mapping:
    # c0 = F2+0.5*(F1+F3)
    # c1 = 0.5*(-F1+3F3)
    # c2 = -F2 + 4F4
    # c3 = -2F3
    # c4 = -4F4
    # c5 = 0
    # see the gromacs manual for more info
    def add_dihedral_RB(self, i, j, k, l, c0, c1, c2, c3, c4, c5):
        self.n_dihedrals += 1
        self.n_params['n_params'][0]['n_dihedrals']= self.n_dihedrals;
        
        self.dihedral_params['dihedrals'].append(OrderedDict([
            ('type', 'RB'),
            ('i', i),
            ('j', j),
            ('k', k),
            ('l', l),
            ('func', 3),
            ('c0', c0),
            ('c1', c1),
            ('c2', c2),
            ('c3', c3),
            ('c4', c4),
            ('c5', c5)
        ]))
    '''Gromacs type 5 (i.e., func =5) Fourier, i.e, OPLS style dihedrals
        note, these do not support phase shifting of phi
        '''
    def add_dihedral_fourier(self, i, j, k, l, c1, c2, c3, c4):
        self.n_dihedrals += 1
        self.n_params['n_params'][0]['n_dihedrals']= self.n_dihedrals;
        
        self.dihedral_params['dihedrals'].append(OrderedDict([
             ('type', 'fourier'),
             ('i', i),
             ('j', j),
             ('k', k),
             ('l', l),
             ('func', 5),
             ('c1', c1),
             ('c2', c2),
             ('c3', c3),
             ('c4', c4)
             ]))
    
    def add_dihedral_periodic(self, i, j, k, l, phi_s, k_phi, multiplicity):
        self.n_dihedrals += 1
        self.n_params['n_params'][0]['n_dihedrals']= self.n_dihedrals;
        
        self.dihedral_params['dihedrals'].append(OrderedDict([
                                                 ('type', 'periodic'),
                                                 ('i', i),
                                                 ('j', j),
                                                 ('k', k),
                                                 ('l', l),
                                                 ('func', 1),
                                                 ('phi_s', phi_s),
                                                 ('k_phi', k_phi),
                                                 ('multiplicity', multiplicity)
                                                 ]))
            
    def add_dihedral_periodic_multiple(self, i, j, k, l, phi_s, k_phi, multiplicity):
        self.n_dihedrals += 1
        self.n_params['n_params'][0]['n_dihedrals']= self.n_dihedrals;
     
        self.dihedral_params['dihedrals'].append(OrderedDict([
                                              ('type', 'periodic multiple'),
                                              ('i', i),
                                              ('j', j),
                                              ('k', k),
                                              ('l', l),
                                              ('func', 9),
                                              ('phi_s', phi_s),
                                              ('k_phi', k_phi),
                                              ('multiplicity', multiplicity)
                                              ]))
    #note, gromacs doesn't differentiate between proper and improper dihedrals in the data file
    #(i.e., both are in the dihedraltypes section)
    #however packages likes lammps put them in different sections
    #so we will store them their own list
    def add_improper_periodic(self, i, j, k, l, phi_s, k_phi, multiplicity):
        self.n_dihedrals += 1
        self.n_params['n_params'][0]['n_impropers']= self.n_impropers;
            
        self.improper_params['impropers'].append(OrderedDict([
             ('type', 'periodic'),
             ('i', i),
             ('j', j),
             ('k', k),
             ('l', l),
             ('func', 4),
             ('phi_s', phi_s),
             ('k_phi', k_phi),
             ('multiplicity', multiplicity)
             ]))
    def print_itp(self, f):
        print >> f, '[ atomtypes ]'
        print >> f, '; name\tbond_type\tatm. nm\tmass\tcharge\tptype\tsigma\tepsilon'
        for atom in self.atom_params['atoms']:
            print >> f,  '%s\t%s\t%d\t%lf\t%lf\t%s\t%lf\t%lf' % (atom['name'], atom['bond_type'], atom['atomic_number'], atom['mass'] , atom['charge'], atom['ptype'], atom['sigma'], atom['epsilon'])

        print >> f, ''
        print >> f, '[ bondtypes ]'
        print >> f, '; i\tj\tfunc\tparameters'
        for bond in self.bond_params['bonds']:
            if (bond['type'] == 'harmonic') or (bond['type'] == 'G96'):
                print >> f, '%s\t%s\t%d\t%lf\t%lf%s%s%s' % (bond['i'], bond['j'], bond['func'], bond['b0'], bond['kb'],' ; ', bond['type'], ' bi j func b0 kb')
            elif bond['type'] == 'morse':
                print >> f, '%s\t%s\t%d\t%lf\t%lf\t%lf\t%s%s%s' % (bond['i'], bond['j'], bond['func'], bond['b0'], bond['D'], bond['beta'], ' ; ', bond['type'], ' i j func b0 D beta')
            elif bond['type'] == 'fene':
                print >> f, '%s\t%s\t%d\t%lf\t%lf\t%s%s%s' % (bond['i'], bond['j'], bond['func'], bond['bm'], bond['kb'], ' ; ', bond['type'], ' i j func bm kb')


        print >> f, ''
        print >> f, '[ angletypes ]'
        print >> f, '; i\tj\tk\tfunc\tth0\t cth'
        for angle in self.angle_params['angles']:
            if (angle['type'] == 'harmonic') or (angle['type'] == 'G96'):
                print >> f, '%s\t%s\t%s\t%d\t%lf\t%lf%s%s%s' % (angle['i'], angle['j'], angle['k'], angle['func'], angle['theta_0'], angle['kb'], ' ; ', angle['type'], ' i j func th0 cth')

        print >> f, ''
        print >> f, '[ dihedraltypes ]'
        print >> f, '; i\tj\tk\tl\tfunc\tcoefficients'
        for dihedral in self.dihedral_params['dihedrals']:
            if dihedral['type'] == 'RB':
                print >> f, '%s\t%s\t%s\t%s\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf' % (dihedral['i'], dihedral['j'], dihedral['k'], dihedral['l'], dihedral['func'], dihedral['c0'], dihedral['c1'],dihedral['c2'],dihedral['c3'],dihedral['c4'],dihedral['c5'])
            elif dihedral['func'] == 'fourier':
                print >> f, '%s\t%s\t%s\t%s\t%d\t%lf\t%lf\t%lf\t%lf' % (dihedral['i'], dihedral['j'], dihedral['k'], dihedral['l'], dihedral['func'], dihedral['c1'], dihedral['c2'],dihedral['c3'],dihedral['c4'])
            elif (dihedral['func'] == 'periodic') or (dihedral['func'] == 'periodic multiple'):
                print >> f, '%s\t%s\t%s\t%s\t%d\t%lf\t%lf\t%lf' % (dihedral['i'], dihedral['j'], dihedral['k'], dihedral['l'], dihedral['func'], dihedral['phi_s'], dihedral['k_phi'],dihedral['multiplicity'])
        for improper in self.improper_params['impropers']:
            if improper['type'] == 'periodic':
                print >> f, '%s\t%s\t%s\t%s\t%d\t%lf\t%lf\t%lf' % (improper['i'], improper['j'], improper['k'], improper['l'], improper['func'], improper['phi_s'], improper['k_phi'],improper['multiplicity'])





    def print_json(self,f):
        forcefield = {}
        forcefield['forcefield'] = []
        forcefield['forcefield'].append(self.n_params)
        forcefield['forcefield'].append(self.atom_params)
        forcefield['forcefield'].append(self.bond_params)
        forcefield['forcefield'].append(self.angle_params)
        forcefield['forcefield'].append(self.dihedral_params)
        forcefield['forcefield'].append(self.improper_params)
        json.dump(forcefield, f)

    def print_json_screen(self):
        forcefield = {}
        forcefield['forcefield'] = []
        forcefield['forcefield'].append(self.n_params)
        forcefield['forcefield'].append(self.atom_params)
        forcefield['forcefield'].append(self.bond_params)
        forcefield['forcefield'].append(self.angle_params)
        forcefield['forcefield'].append(self.dihedral_params)
        forcefield['forcefield'].append(self.improper_params)
        print json.dumps(forcefield, indent = 2)
    
    def load_json(self, f):
        forcefield = {}
        forcefield['forcefield'] = []
        forcefield = json.load(f)
        
        #first load in the number of each type of parameter
        #this obviously isn't explicitly needed, but is still useful to know and have available to us
        for i in range(0, len(forcefield['forcefield'])):
            if 'n_params' in forcefield['forcefield'][i].keys():
                for params in forcefield['forcefield'][i]['n_params']:
                    self.n_atoms = params['n_atoms']
                    self.n_bonds = params['n_bonds']
                    self.n_angles = params['n_angles']
                    self.n_dihedrals = params['n_dihedrals']
                    self.n_impropers = params['n_impropers']

        for i in range(0, len(forcefield['forcefield'])):
            if self.n_atoms > 0:  #not really necessary because of the if statement below, but this avoids doing extra work
                if 'atoms' in forcefield['forcefield'][i].keys():
                    for atom in forcefield['forcefield'][i]['atoms']:
                        self.atom_params['atoms'].append(atom)
            if self.n_bonds > 0:
                if 'bonds' in forcefield['forcefield'][i].keys():
                    for bond in forcefield['forcefield'][i]['bonds']:
                        self.bond_params['bonds'].append(bond)
            
            if self.n_angles > 0:
                if 'angles' in forcefield['forcefield'][i].keys():
                    for angle in forcefield['forcefield'][i]['angles']:
                        self.angle_params['angles'].append(angle)
            if self.n_dihedrals > 0:
                if 'dihedrals' in forcefield['forcefield'][i].keys():
                    for dihedral in forcefield['forcefield'][i]['dihedrals']:
                        self.dihedral_params['dihedrals'].append(dihedral)
            if self.n_impropers > 0:
                if 'impropers' in forcefield['forcefield'][i].keys():
                    for improper in forcefield['forcefield'][i]['impropers']:
                        self.improper_params['impropers'].append(improper)

def main():
    
    out_file_itp = open('parameters.itp','w')
    out_file_json = open('parameters.json','w')
    params = Parameters()

    #note gromacs itp format expects nm for sigma and kJ/mol for epsilon
    params.add_atom(name='opls_135', bond_type='CT', atomic_number=6, mass=12.01100, charge=-0.18, ptype='A', sigma=0.35, epsilon=0.276144)
    params.add_atom(name='opls_136', bond_type='CT', atomic_number=6, mass=12.01100, charge=-0.12, ptype='A', sigma=0.35, epsilon=0.276144)
    params.add_atom(name='opls_140', bond_type='HC', atomic_number=1, mass=1.00800, charge=0.06, ptype='A', sigma=0.25, epsilon=0.125520)

    #add a harmonic bond. Note gromacs has a 1/2 prefactor, lammps does not and the original OPLS definition doe snot.
    params.add_bond_harmonic(i='CT', j='CT', b0=0.15290, kb=224262.4)

    #note GROMACS has a factor of 1/2 in front of the harmonic angle term: lammps does not, OPLS original definition does not,
    #i.e., energy parameter kb is twice as large as in the original opls paper
    params.add_angle_harmonic(i='CT', j='CT', k='CT', theta_0=112.700, kb=488.273)


    #add a dihedral.  These are OPLS style fourier dihedrals, recast into RB form
    params.add_dihedral_RB(i='CT', j='CT', k='CT', l='CT', c0=2.92880, c1=-1.46440, c2=0.20920, c3=-1.67360, c4=0.00000, c5=0.00000)
    
    #print the itp and json files
    params.print_itp(out_file_itp)
    params.print_json(out_file_json)

    params.print_json_screen()

    out_file_json.close()
    
    
    #make sure we can properly read the json file
    in_file_json = open('parameters.json','r')

    nbParams2 = Parameters()

    nbParams2.load_json(in_file_json)

    out_file_itp2 = open('parameters2.itp','w')
    
    #this itp file should match the one from earlier
    nbParams2.print_itp(out_file_itp2)


if __name__ == "__main__":
    main()
