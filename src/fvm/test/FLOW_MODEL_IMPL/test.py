# redirect stdout to a disk file
import sys
saveout = sys.stdout
outfile = open('output.txt', 'w')
sys.stdout = outfile

# output stuff
print 'hello world'

# restore stdout
outfile.flush()
outfile.close()
sys.stdout = saveout
