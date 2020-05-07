import swc_intersect_detect
from swc_intersect_detect import morph_io

# Convert nrn to swc file
# Not necessary if the acquired dat is in .swc format
morph_io.nrn2swc("Purkinje19b972-1.nrn", "Purkinje.swc")

# perform intersection detection
swc_intersect_detect.run("Purkinje.swc")