import os
from matplotlib.backends.backend_pdf import PdfPages

from . import Step01
from . import Step02
from . import Step03
from . import Step04
from . import Step05
from . import Step06
from . import Step07
from . import Step08
from . import Step09
from . import Step10
from . import Step11
from . import Step12

output = os.path.join(os.path.dirname(__file__), "Gallery.pdf")

with PdfPages(output) as pdf:
    for index, step in enumerate([
            Step01, Step02, Step03, Step04, 
            Step05, Step06, Step07, Step08, 
            Step09, Step10, Step11, Step12]):
    
        step.main()
        pdf.savefig()
