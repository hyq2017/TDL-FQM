import os
import numpy as np

def input_coord(coord, atoms):
    a = np.loadtxt(coord)
    b = np.reshape(a,(1, atoms * 3))
    return b

# read the type and index of atoms
def type_index():
    f = open('type.raw','r')
    for line in f:
        types = np.array([int(x) for x in line.split()])
    f.close()
    index_C = np.where(types == 0)
    index_H = np.where(types == 1)
    index_O = np.where(types == 2)
    index_N = np.where(types == 3)
    index_S = np.where(types == 4)
    index = np.concatenate((index_C,index_H,index_O,index_N,index_S),axis=1)
    return index

def main():
    #os.system("cp ../package/ALA/type.raw .")
    #os.system("cp ../input.ALA .")
    atoms = 29
    i = 0
    energy_bias = 0.0001
    f_box = open('box.raw', 'a+')
    input_file = list()
    a = np.zeros((1, atoms * 3))
    for line in open("input.HIE"):
        fragment = line.strip('\n')
        os.system("cp ../" + fragment + " ." )
        input_file.append(line.strip('\n'))
        b = input_coord(fragment, atoms)
        a = np.concatenate((a, b), axis=0)
        i = i + 1
        f_box.write('100.0 0.0 0.0 0.0 100.0 0.0 0.0 0.0 100.0' + '\n')
    np.savetxt('coord.raw', a[1:])
    f_box.close()

    os.system('../package/bin/raw_to_set.sh ' + str(i))
    os.system('../package/bin/dp_test -m ../package/model/HIE.pb -s . -S set -n ' + str(i) + ' -d detail')
    os.system('rm box.raw')

    # read energy
    f = open('detail.e.out', 'r')
    next(f)
    energy = list()
    for line in f:
        energy.append(np.array([float(x) for x in line.split()])[1])
    f.close()
    
    # write energy
    for i in range(len(input_file)):
        f_e = open(input_file[i][0:len(input_file[i])-3] + 'log', 'a+')
        f_e.write('energy: ' + '\n' + str(energy[i] / 27.21138602 - energy_bias))
        f_e.write('\n')
        f_e.close()

    # read disorder force
    f = open('detail.f.out', 'r')
    next(f)
    force_disorder = np.zeros((atoms * 100, 3))
    num = 0
    for line in f:
        force_disorder[num,:] = np.array([float(x) for x in line.split()])[3:6]
        num=num+1
    f.close()

    # read the type and index of atoms
    index = type_index()
  

    # get order force
    count = int(num / atoms)
    force_sub  = np.zeros((count, atoms, 3))
    for i in range(count):
        force_sub[i][0:atoms, :] = force_disorder[(i*atoms):((i+1)*atoms), :] / 27.21138602 * 0.529177249
    
    force = np.zeros((count, atoms, 3))
    for i in range(len(input_file)):
        count = 0
        for j in index[count]:
            force[i][j,:] = force_sub[i][count,:]
            count = count + 1
        f_f = open(input_file[i][0:len(input_file[i])-3] + 'log', 'a+')
        f_f.write('force:' + '\n')
        f_f.write(str(force[i]).replace('[', '').replace(']', '') + '\n')
        f_f.close()

    #os.system("rm detail.*")
    #os.system("rm coord.raw")
    #os.system("rm -rf set.000")

if __name__ == "__main__":
    main()
