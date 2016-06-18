import argparse


def fixply(filename, outputname):
    num_lines = 0
    with open(filename) as pointsfile:
        line = pointsfile.readline()
        while(line != ""):
            num_lines += 1
            line = pointsfile.readline()
    with open(outputname, "w") as fixedfile:
        fixedfile.write(
'''ply
format ascii 1.0
comment VCGLIB generated
element vertex {0}
property float x
property float y
property float z
property uchar red
property uchar green
property uchar blue
end_header
'''.format(num_lines))
        with open(filename) as pointsfile:
            line = pointsfile.readline()
            while(line != ""):
                fixedfile.write(line)
                line = pointsfile.readline()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    This script reads a list of 3D points (x, y) meant to be in a Ply format and adds a Ply header
    ''')
    parser.add_argument('input_file', help='input list of points per line (format: ply)')
    parser.add_argument('output_file', help='name of the output file')
    args = parser.parse_args()
    fixply(args.input_file, args.output_file)
