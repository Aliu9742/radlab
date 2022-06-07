"""
This is a list of functions that was created
for the BF-II/membrane binding project, but
they are intended to be reusable
for future Radlab members regardless of the project

A lot of the functions here have
only gone through limited testing,
so I definitely recommend printing out what
they return to make sure they match your expectations.

It is my hope that future lab members will 
continue to use and refine the functions here.  

For questions please reach out to aliu9@wellesley.edu
or 425-559-4470
"""

import os
import copy

###############################################################################################


class Atom():
    def __init__(self, line):
        """
        pdb, crd, or gro file lines are accepted when creating an Atom object
        both pdb and crd file use Angstroms for position but gro file uses nm
        when position of gro files are read, it is automatically converted to units of Angstroms
        """

        line = line.strip('\n')
        if line[:4] == 'ATOM':
            # reading pdb file
            attr_dict = self.pdb_to_dict(line)

        elif line[:6].strip().isdigit() and line[6:11].strip().isdigit:
            # reading crd file assumes first 2 columns are numbers
            attr_dict = self.crd_to_dict(line)

        elif line[:5].strip().isdigit() and (not line[5:11].strip().isdigit()):
            # reading gro file assumes the first column is a number and the second column is not a number
            attr_dict = self.gro_to_dict(line)

        # if attributes are missing, insert default
        defaults = {'res_num': '',
                    'chain_id': '',
                    'alternate_location': '',
                    'res_code': '',
                    'element': '',
                    'occupancy': '1.00',
                    'seg_identifier': '',
                    'chain': '',
                    'charge': '',
                    'is_amino_acid': is_amino_acid(attr_dict['res_name'])
                    }

        for attr, val in defaults.items():
            if attr not in attr_dict:
                attr_dict[attr] = val

        # convert the attribute dictionary to the attribute of Atom object
        for attr in attr_dict:
            setattr(self, attr, attr_dict[attr])

    def pdb_to_dict(self, line):
        """
        create a dictionary by reading in line from pdb file
        with keys:
        atom_num, atom_name, alternate_location, res_name
        chain, res_num, res_code, xpos, ypos, zpos (in Angstroms)
        occupancy, temp_factor, seg_identifier, element
        all values are strings

        cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html 
        for structure of pdb file
        """

        attr_dict = {'num': line[6:11].strip(),

                     # .strip() removes extra spaces before and after a string
                     # name doesn't get stripped to preserve the spacing of columns when converting back to pdb line
                     'name': line[12:16],
                     'alternate_location': line[16],

                     # difference from website
                     # lipids POPG/POPE can be 4 letters instead of 3
                     'res_name': line[17:21].strip(),

                     'chain': line[21],
                     'res_num': line[22:26].strip(),
                     'res_code': line[26],
                     'xpos': line[30:38].strip(),
                     'ypos': line[38:46].strip(),
                     'zpos': line[46:54].strip(),
                     'occupancy': line[54:60].strip(),

                     # difference from website
                     # radlab uses temp factor for charge
                     'charge': line[60:66].strip(),
                     'seg_identifier': line[72:76].strip(),
                     'element': line[76:78].strip()
                     }
        return attr_dict

    def crd_to_dict(self, line):
        """
        Take in a line from crd file return a dictionary
        with keys:
        atom_num, res_num, res_name, atom_name
        xpos, ypos, zpos (in Angstroms), chain, chain_id, charge 

        all values are strings
        """

        attr_dict = {
            'num': line[:6].strip(),
            'res_num': line[6:11].strip(),
            'res_name': line[11:16].strip(),
            'name': line[16:20].strip(),
            'xpos': line[21:30].strip(),
            'ypos': line[32:40].strip(),
            'zpos': line[41:50].strip(),
            'chain': line[51:52].strip(),
            'chain_id': line[55:61].strip(),
            'charge': line[61:].strip()
        }
        return attr_dict

    def gro_to_dict(self, line):
        """
        Take in an atom line from gro file
        return a dictionary 
        res, atom, atom_num, xpos, ypos, zpos

        positions are converted from nm to Angstroms

        for format of gro file see
        https://manual.gromacs.org/archive/5.0.3/online/gro.html
        """
        def to_angstrom(num): return str(float(num.strip())*10)

        attr_dict = {
            'res_num': line[:5].strip(),
            'res_name': line[5:11].strip(),
            'name': line[11:16].strip(),
            'num': line[16:21].strip(),
            'xpos': to_angstrom(line[21:29]),
            'ypos': to_angstrom(line[29:37]),
            'zpos': to_angstrom(line[37:45]),
        }

        return attr_dict

    def pdb_line(self):
        """
        format atom's attributes into pdb line return as string
        referenced https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
        for structure of pdb

        the atom names might be misaligned but the pdb file should still load in vmd
        """
        floats = {'xpos': self.xpos,
                  'ypos': self.ypos,
                  'zpos': self.zpos,
                  'occupancy': self.occupancy,
                  'charge': self.charge}
        floats_formatted = {}
        def pos_format(a): return '{:>8.3f}'.format(float(a))
        def occu_charge_format(a): return '{:>6.2f}'.format(float(a))

        for attr, val in floats.items():
            if 'pos' in attr:
                try:
                    floats_formatted[attr] = pos_format(val)
                except:
                    floats_formatted[attr] = ''
            else:
                try:
                    floats_formatted[attr] = occu_charge_format(val)
                except:
                    floats_formatted[attr] = ''

        return 'ATOM  {:>5} {:4}{:1}{:<4}{:1}{:>4}{:1}   {:>8}{:>8}{:>8}{:>6}{:>6}      {:<4}{:>2}\n'.format(
            self.num,
            self.name,
            self.alternate_location,
            self.res_name,
            self.chain,
            self.res_num,
            self.res_code,
            floats_formatted['xpos'],
            floats_formatted['ypos'],
            floats_formatted['zpos'],
            floats_formatted['occupancy'],
            floats_formatted['charge'],
            self.seg_identifier,
            self.element)

    def crd_line(self):
        """
        format atom's attributes to a crd line
        """
        floats = {'xpos': self.xpos,
                  'ypos': self.ypos,
                  'zpos': self.zpos,
                  'charge': self.charge}

        floats_2 = {}

        for attr, val in floats.items():
            try:
                floats_2[attr] = '{:9.5f}'.format(float(val))
            except:
                floats_2[attr] = ''

        return '{:>5}{:>5}{:>4}  {:<4} {:10}{:10}{:10}{:<4}{:>4}  {:10}\n' .format(
            self.num,
            self.res_num,
            self.res_name,
            self.name.strip(),
            floats_2['xpos'],
            floats_2['ypos'],
            floats_2['zpos'],
            self.chain,
            self.chain_id,
            floats_2['charge'])

    def gro_line(self):
        """
        format atom's attributes to a gro line
        https://manual.gromacs.org/archive/5.0.3/online/gro.html
        for format of gro file
        the last three columns with velocities are omitted
        """

        def to_nm(num): return float(num.strip())/10

        x = '{:8.3f}'.format(to_nm(self.xpos))
        y = '{:8.3f}'.format(to_nm(self.ypos))
        z = '{:8.3f}'.format(to_nm(self.zpos))

        return '{:>5}{:<5}{:>5}{:>5}{}{}{}\n'.format(
            self.res_num,
            self.res_name,
            self.name,
            self.num,
            x,
            y,
            z
        )

###############################################################################################

def crd_line_contains_atom(line):
    """
    If crd line is header line, 
    return false, 
    if it contains information about atom, 
    return true
    """

    if line[0] == '*' or line.strip().isdigit():
        return False
    else:
        return True

###############################################################################################


def get_pdb_box_size(line):
    """
    take in pdb line, return box dimension as tuple of floats (x,y,z)
    """
    line_split = line.strip().split()
    return float(line_split[1]), float(line_split[2]), float(line_split[3])

###############################################################################################

def change_atom_name(atom_name):
    """
    take in atom_name as str and moves the first character of of the atom_name to
    the last position eg. given 1HG return HG1, 12HG return HG12 

    This function is needed because there's a bug in our LPBE solver that can't have
    atom_names starting with a number, so during snapshot preparation, atom_names are changed
    """
    if atom_name.isdigit():
        raise ValueError()

    if not atom_name[0].isdigit():
        return atom_name
    else:
        return change_atom_name(atom_name[1:] + atom_name[0])

####################################################################################

def is_amino_acid(res_type):
    """
    Given 3 letter string res_type (capitalization ignored) 
    return true if residue is amino acid and 
    false is res_type is not amino acid (sol, lipid)
    Raise exception if res_type is anything other than letters
    """

    if not isinstance(res_type, str):
        raise TypeError(
            'is_amino_acid takes str argument, got {}'.format(res_type))

    amino_acids = ['arg', 'his', 'lys', 'asp', 'glu',
                   'ser', 'thr', 'asn', 'gln', 'cys',
                   'sec', 'gly', 'pro', 'ala', 'val',
                   'ile', 'leu', 'met', 'phe', 'tyr', 'trp']

    if res_type.lower() in amino_acids:
        return True
    else:
        return False

#############################################################################

def only_alphabet(string):
    """
    Given a string, return only the characters that are alphabets
    This function was created because the res column in a gro file
    has both res num and res type eg. 12ARG
    This function could be used to extract res type eg. ARG
    """
    string_list = [char for char in string if char.isalpha()]
    return ''.join(string_list)

#############################################################################

def only_digit(string):
    """
    Given a string, return only the characters that are numbers
    eg only_digit(file13.crd) will return 13
    only_digit(1file3.crd) will return 13 too
    """
    string_list = [char for char in string if char.isdigit()]
    return ''.join(string_list)

##############################################################################

def read_itp_atom_lines(line):
    """
    given line from [ atom ] section of itp file
    return None if line starts with ;
    return dictionary with keys atom_type, res_id, res_type, atom_name, charge 
    """
    if not isinstance(line, str):
        raise TypeError(
            'read_itp_atm_lines takes str argument, got {}'.format(type(line)))

    if line[0] != ';':
        line_split = line.split()
        line_dict = {
            'atom_type': line_split[1],
            'res_id': int(line_split[2]),
            'res_type': line_split[3],
            'atom_name': line_split[4],
            'charge': float(line_split[6])
        }
        return line_dict

##################################################################################################

def itp_atom_lines(f):
    """
    return [ atoms ] section of itp file as a list of strings
    by looking for lines that occur between [ atoms ] and [ bonds ] header
    to create a list called itp_lines
    
    using itp_lines list
    returns a dictionary with atom_name : charge as key val pairs 
    this function was written for adding charges into crd files
    """
    if not isinstance(f, str):
        raise TypeError(
            'itp_file_lines takes str argument, got {}'.format(type(f)))
    with open(f, 'r') as itp_file:
        itp_lines = [line.strip() for line in itp_file]

        for index, line in enumerate(itp_lines):
            if line == '[ atoms ]':
                top_line_index = index
            if line == '[ bonds ]':
                bottom_line_index = index

    itp_lines = [line for index, line in enumerate(itp_lines) if index > top_line_index+1 and index < bottom_line_index-1]
    itp = (read_itp_atom_lines(i) for i in itp_lines if read_itp_atom_lines(i))
    return {i['atom_name']: i['charge'] for i in itp} 

############################################################################################

def itp_total_charge(atom_charges):
    """
    given a dictionary with atom_name and charges as key value pairs,
    (output of itp_atom_lines)
    return net charge of molecule
    """
    sum = 0
    for i in atom_charges.values():
        sum += i
    return sum   

###########################################################################################

def return_itp_molecule_type(itp):
    """
    given path of itp file, open and read it and
    return molecule type eg (POPE, POPG, Protein)
    """
    if not isinstance(itp, str):
        raise TypeError(
            'return_itp_molecule_type takes str argument, got {}'.format(type(itp)))
    with open(itp, 'r') as in_file:
        counter = 0
        found_moleculetype = False
        for line in in_file:
            if line.strip() == '[ moleculetype ]':
                found_moleculetype = True
            if found_moleculetype:
                counter += 1
            if counter == 3:
                mole_type = line.split()[0].strip()
                return mole_type

############################################################################################

def create_gams_table(header, set_i, set_j, matrix, col_spacing=15):
    """
    take in header as string: 'd(i,j) distance in thousands of miles'
    take in set_i as array of strings: ['seattle', 'san-diego']
    take in set_j as array of strings: ['new-york', 'chicago', 'topeka']
    take in matrix as 2d array of strings : [['2.5', '1.7', '1.8'], ['2.5', '1.8', '1.4']]
    col_spacing optional argument, default 15, adjust bigger number for wider spacing, smaller number for tighter spacing 
    return str:
    Table
    d(i,j) distance in thousands of miles
            new-york chicago topeka
    seattle 2.5       1.7     1.8
    san-diego 2.5     1.8     1.4

    assumes set_i is the rows, set_j is col
    Note this function doesn't work if
    numpy array is fed in instead of python array
    """
    set1 = copy.deepcopy(set_i)
    set2 = copy.deepcopy(set_j)
    m = copy.deepcopy(matrix)

    rows = ['']*len(m)

    for i, row in enumerate(m):
        m[i].insert(0, set1[i])

        for j, col in enumerate(row):
            row[j] = '{:{}}'.format(col, col_spacing)

        rows[i] = ''.join(m[i])

    for index, j in enumerate(set2):
        set2[index] = '{:{}}'.format(j, col_spacing)

    set2.insert(0, ' '*col_spacing)
    table_col_header = "".join(set2)

    table_rows = '\n'.join(rows)

    return "Table\n{}\n{}\n{} ;\n".format(header, table_col_header, table_rows)

################################################################################################

def create_table(matrix, col_spacing=15):
    """
    take in matrix output a table as string
    col_spacing optional argument, default 15, adjust bigger number for wider spacing, smaller number for tighter spacing
    this function is different from create_gams_table in that
    it doesn't take in headers
    """
    rows = [0] * len(matrix)
    for index, row in enumerate(matrix):
        if isinstance(row, list):
            l = ['{:{}}'.format(col, col_spacing) for col in row]
            rows[index] = ''.join(l)
        else:
            rows[index] = row

    table_joined = '\n'.join(rows)

    return table_joined

###############################################################################

def create_gams_list(header, set_i, vector):
    """
    take in header as string: 'a(i) capacity of plant i in cases'
    take in set_ i as an array of strings: ['seattle', 'san-diego']
    take in vector as an array of strings: ['350', '600']

    return str:
    Parameters
    a(i) capacity of plant i in cases
    /	seattle	350
            san-diego	600	/;
    """
    list_rows = [0] * len(set_i)

    # print(set_i)
    # print(vector)
    for i, ele in enumerate(set_i):
        list_rows[i] = '\t'.join((ele, vector[i]))

    rows = '\n\t\t'.join(list_rows)
    #print('Parameter\n\t{}\n\t\t/{}/ ;\n'.format(header, rows))
    return 'Parameter\n\t{}\n\t\t/{}/ ;\n'.format(header, rows)

###########################################################################

def create_gams_scalar(name_description, val):
    """
    take in name_description as string: 'f freight in dollars per case per thousand miles
    take in val as string or float or int: 90

    return str:
    Scalar f frieght in dollars per case per thousand miles /90/;
    """

    return '\nScalar {} /{}/ ;\n'.format(name_description, val)

##############################################################################

def create_gams_sets(sets):
    """
    take in sets as dictionary, the key is name and description of the set
    the value is the members of the set

    {'i canning plants' : ['seattle', 'san-diego'],
     'j markets' : ['new-york', 'chicago', 'topeka'],
    }
    return str:
    Sets
            i canning plants	/seattle, san-diego/
            j markets	/new-york, chicago, topeka/;
    """
    row = [''] * len(sets)

    for index, (key, val) in enumerate((sets.items())):
        set_members = ['/'] * (len(val)+2)

        for i, member in enumerate(val, start=1):
            if len(val) == i:
                set_members[i] = '{} '.format(member)
            else:
                set_members[i] = '{}, '.format(member)

        set_members_str = ' '.join(set_members)

        row[index] = '\t'.join([key, set_members_str])

    rows = '\n\t'.join(row)

    #print('\nSets\n\t{} ;\n'.format(rows))
    return '\nSets\n\t{} ;\n'.format(rows)

##################################################################################

def extract_from_gams(gams_file, var, cons=True, float_out=False):
    # This function is incomplete
    """
    given path of gams_file and name of variable to look for in gams file (output of gams with extension lst)
    if cons is True (default value) return a float eg. deltaG
    else return an array eg. optimal qs

    float_out set to False by default, will return as strs or list of strs

    This function makes the assumption that the output lines is unbroken:
    ---- VAR qs  charges on side chain

          LOWER          LEVEL          UPPER         MARGINAL

q1        -1.0000        -1.0000         1.0000         2.0401      
q2        -1.0000         0.3008         1.0000          .          
q3        -1.0000         1.0000         1.0000        -0.0255      
q4        -1.0000         0.3353         1.0000         EPS         
q5        -1.0000        -0.6835         1.0000   -1.254636E-9      
q6        -1.0000         1.0000         1.0000        -3.0361      
q7        -1.0000        -0.6557         1.0000         EPS         
q8        -1.0000         0.7477         1.0000         EPS
    instead of:
    ---- VAR qs  charges on side chain

          LOWER          LEVEL          UPPER         MARGINAL

q1        -1.0000        -1.0000         1.0000         2.0401      
q2        -1.0000         0.3008         1.0000          .          
q3        -1.0000         1.0000         1.0000        -0.0255      
q4        -1.0000         0.3353         1.0000         EPS         
q5        -1.0000        -0.6835         1.0000   -1.254636E-9
----------------page break     
q6        -1.0000         1.0000         1.0000        -3.0361      
q7        -1.0000        -0.6557         1.0000         EPS         
q8        -1.0000         0.7477         1.0000         EPS         


    """
    with open(gams_file, 'r') as f:

        if cons:
            for line in f:
                if '---- VAR {}'.format(var) in line:
                    if float_out:
                        return float(line.split()[4])
                    else:
                        return line.split()[4]

        else:
            lines = []
            append = False
            count = 0

            for line in f:
                if append:
                    lines.append(line)
                if '---- VAR {}'.format(var) in line:
                    append = True

                if append and line == '\n':
                    count += 1
                    if count == 3:
                        break

            if float_out:
                return [float(i.split()[2]) for i in lines if i != '\n' and not ('LEVEL' in i)]
            else:
                return [i.split()[2] for i in lines if i != '\n' and not ('LEVEL' in i)]

##############################################################################################################

def generate_run_file_output(atoms_charged, name=0):
    """
    given the atoms being charged, and the name of the grp return a string of 3 lines:
    generate_run_file_output('0', '0')
    will return:
    mark=output
    atoms_charged=(atomno eq "0")
    name="0"

    name is an optional argument, if no argument given, name is assumed to be the same as atoms_charged
    """

    if not(name):
        name = atoms_charged

    if str(atoms_charged).isdigit():
        return 'mark=output\natoms_charged=(atomno eq "{}")\nname="{}"\n\n'.format(atoms_charged, name)
    else:
        return 'mark=output\natoms_charged=(segid eq "{}")\nname="{}"\n\n'.format(atoms_charged, name)

#####################################################################################

def generate_run_file_bound_unbound_state(mark, state, atoms_charged, atoms_shape, atoms_center):
    # This function could be made more user friendly so the user doesn't have to put in the () or the "" marks
    """
    given 4 strs return a group of 4 lines in the run file,
    user need to provide parenthesis and quotation marks 
    """
    return 'mark={}\nstate="{}"\natoms_charged={}\natoms_shape={}\natoms_center={}\n\n'.format(mark, state, atoms_charged, atoms_shape, atoms_center)

#########################################################################################

def create_run_file_partial_charge_opt(s_atoms, current_atom):
    """
    s_atoms is a list of side chain atoms
    return entirety of run file
    """
    start = 'mark=start\n\n'
    end = 'mark=end'

    side_chain = ''
    for s in s_atoms:
        side_chain += generate_run_file_output(s)

    Lstqt = generate_run_file_output('T', 'Lstqt')
    Cqr = generate_run_file_output('R', 'Cqr')

    bound = generate_run_file_bound_unbound_state('final', 'bound_L', '(atomno eq "{}")'.format(
        current_atom), 'all', '(atomno eq "{}")'.format(current_atom))
    unbound = generate_run_file_bound_unbound_state('reference', 'unbound_L', '(atomno eq "{}")'.format(
        current_atom), '(segid eq "S" || segid eq "T")', '(atomno eq "{}")'.format(current_atom))

    return start + side_chain + Lstqt + Cqr + bound + unbound + end

#############################################################################################

def read_diff_table_w_4_vals(diff_table_path):
    """
    Given path of diff table (for deltaG calculation with 4 values not for 
    charge optimization), open and read the file
    return a dictionary with keys RDP, LDP, and INT
    vals are numbers as floats
    """

    with open(diff_table_path, 'r') as f:
        lines = [i.strip().split() for i in f.readlines()]

    deltaG_comps = {'LDP': float(lines[1][1]),
                    'INT': (float(lines[1][4])+float(lines[2][1])) / 2,
                    'RDP': float(lines[2][4])
                    }
    return deltaG_comps

#################################################################################################

def extract_matrix_vec_cons(res_dir, vec1='Cqr', vec2='Lstqt', cons1='Ltt', cons2='Cqrt', cons_dir='constants'):
    """
    given the path to res_dir as string, assemble Lmatirx, Cqr vector, and Lstqt vector
    and return as dictionary {Lmatrix : [], Cqr : [], Lstqt : []...}
    numbers in arrays are of type string
    Takes in optional parameters for name of output grp, default is Cqr and Lstqt vectors
    assumes atom dirs are numbered, and looks for atoms in difference.table file as numbers 
    """

    def str_to_int(string):
        return int(string)

    #find all the numbered directories in given res_dir and sort them from smallest to larges
    atom_dir = [i for i in os.listdir(res_dir) if os.path.isdir(
        os.path.join(res_dir, i)) and i.isdigit()]
    atom_dir.sort(key=str_to_int)

    #create dictionary for filling in values in the rest of the function
    ele = {'Lmatrix': [[0 for x in range(len(atom_dir))] for y in range(len(atom_dir))],
           vec1: [0] * len(atom_dir),
           vec2: [0] * len(atom_dir),
           cons1: '0.00',
           cons2: '0.00'
           }

    for Lindex, atom in enumerate(atom_dir):
        #reading difference.table to get matrix and vector values
        diff_table = os.path.join(res_dir, atom, 'difference.table')

        with open(diff_table, 'r') as f:
            table = [i.strip().split()
                     for index, i in enumerate(f.readlines()) if index != 0]

            j = 0
            for line in table:
                if line[0].isdigit():
                    if line[0] == atom:
                        ele['Lmatrix'][Lindex][j] = line[1]
                    else:
                        ele['Lmatrix'][Lindex][j] = str(float(line[1])/2)

                    j += 1

                elif vec1 == line[0]:
                    ele[vec1][Lindex] = line[1]
                elif vec2 == line[0]:
                    ele[vec2][Lindex] = line[1]

    #getting Ltt and Cqrt
    diff_table = os.path.join(res_dir, cons_dir, 'difference.table')
    with open(diff_table, 'r') as f:
        table = [i.strip().split()
                 for index, i in enumerate(f.readlines()) if index != 0]
        for line in table:
            if line[0] == cons1:
                ele[cons1] = line[1]
            if line[0] == cons2:
                ele[cons2] = line[1]

    return ele

###################################################################################