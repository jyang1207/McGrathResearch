
# use tkinter for GUI elements
import types
from optparse import OptionParser
import os
import simUtility
import sumSimUtility
import plotSimUtility

try:
    import tkinter as tk # python 3
    tkf = tk.filedialog
    tkm = tk.messagebox
except ImportError:
    import Tkinter as tk # python 2
    import tkFileDialog as tkf
    import tkMessageBox as tkm
    
class ModelerDashboard:
    '''
    class to build and manage the display
    '''
    
    def __init__(self, width, height, verbose=True):
        '''
        Constructor for the display application
        '''
        # make initial fields
        self.verbose = verbose
        self.images = None
        self.runDirectory = None
        self.results = None
        self.summary = None
        
        # create a tk object, which is the root window
        self.root = tk.Tk()

        # set up the geometry for the window
        self.root.geometry( "%dx%d+50+30" % (width, height) )
        
        # set the title of the window
        self.root.title("Modeling Dashboard")
        # set the maximum size of the window for resizing
        self.root.maxsize( 1600, 900 )
        
        # setup the menus
        self.buildMenus()

        # build the controls
        self.buildControls()
        
        # build the status
        self.buildStatus()

        # bring the window to the front
        self.root.lift()

        # - do idle events here to get actual canvas size
        self.root.update_idletasks()

        # set up the key bindings
        self.setBindings()
    
    def buildControls(self):
        '''
        build the frame and controls for application
        '''
        if self.verbose: print("building the control frame")
        # make the control frames
        headers = ["model", "summarize", "display", "plot"]
        headers = [header.capitalize()+"\n" for header in headers]
        frames = [tk.Frame(self.root) for _ in range(len(headers))]
        self.modelFrame, self.sumFrame, self.displayFrame, self.plotFrame = frames
        
        # pack all frames with separators and headers
        for header, frame in zip(headers, frames):
            frame.pack(side=tk.LEFT, padx=2, pady=2, fill=tk.Y)
            tk.Frame( self.root, width=2, bd=1, relief=tk.SUNKEN
                      ).pack( side=tk.LEFT, padx = 2, pady = 2, fill=tk.Y)
            tk.Label( frame, text=header
                      ).pack(side=tk.TOP)
                      
        bwidth = 20
        # make the model frame controls
        tk.Button( self.modelFrame, text="Select Images", 
                   command=self.selectImages, width=bwidth
                   ).pack(side=tk.TOP)
        tk.Button( self.modelFrame, text="Select Run Directory", 
                   command=self.selectRunModelDir, width=bwidth
                   ).pack(side=tk.TOP)
        self.modelBulgeOpt = tk.IntVar()
        tk.Checkbutton( self.modelFrame, text="include bulge component",
                        variable=self.modelBulgeOpt
                        ).pack(side=tk.TOP)
        self.modelGalfitOffOpt = tk.IntVar()
        tk.Checkbutton( self.modelFrame, text="dont run GALFIT",
                        variable=self.modelGalfitOffOpt
                        ).pack(side=tk.TOP)
        self.modelParallelOpt = tk.IntVar()
        tk.Checkbutton( self.modelFrame, text="run in parallel",
                        variable=self.modelParallelOpt
                        ).pack(side=tk.TOP)
        self.modelRealOpt = tk.IntVar()
        tk.Checkbutton( self.modelFrame, text="images degraded",
                        variable=self.modelRealOpt
                        ).pack(side=tk.TOP)
        tk.Label( self.modelFrame, text="MPZ" ).pack(side=tk.TOP)
        self.modelMPZEntry = tk.Entry( self.modelFrame, width=bwidth )
        self.modelMPZEntry.pack(side=tk.TOP)
        tk.Label( self.modelFrame, text="Plate Scale" ).pack(side=tk.TOP)
        self.modelPlateEntry = tk.Entry( self.modelFrame, width=bwidth )
        self.modelPlateEntry.pack(side=tk.TOP)
        tk.Button( self.modelFrame, text="Run Modeling", 
                   command=self.runModel, width=bwidth
                   ).pack(side=tk.TOP) 
        
        # make the summary frame controls
        tk.Button( self.sumFrame, text="Select Results", 
                   command=self.selectResults, width=bwidth
                   ).pack(side=tk.TOP)
        self.sumBulgeOpt = tk.IntVar()
        tk.Checkbutton( self.sumFrame, text="bulge components included",
                        variable=self.sumBulgeOpt
                        ).pack(side=tk.TOP)
        self.sumRealOpt = tk.IntVar()
        tk.Checkbutton( self.sumFrame, text="images degraded",
                        variable=self.sumRealOpt
                        ).pack(side=tk.TOP)
        tk.Label( self.sumFrame, text="Delimiter" ).pack(side=tk.TOP)
        self.sumDelimEntry = tk.Entry( self.sumFrame, width=bwidth
                                       ).pack(side=tk.TOP)
        tk.Label( self.sumFrame, text="Output Filename" ).pack(side=tk.TOP)
        self.sumOutputEntry = tk.Entry( self.sumFrame, width=bwidth
                                       ).pack(side=tk.TOP)
        tk.Button( self.sumFrame, text="Run Summarizing", 
                   command=self.runSummary, width=bwidth
                   ).pack(side=tk.TOP)
        
        # make the display frame controls
        tk.Button( self.displayFrame, text="Select Summary", 
                   command=self.selectSummary, width=bwidth
                   ).pack(side=tk.TOP)
        tk.Button( self.displayFrame, text="Run Displaying", 
                   command=self.runDisplay, width=bwidth
                   ).pack(side=tk.TOP)
        
        # make the plot frame controls
        tk.Button( self.plotFrame, text="Select Summary", 
                   command=self.selectSummary, width=bwidth
                   ).pack(side=tk.TOP)
        tk.Button( self.plotFrame, text="Run Plotting", 
                   command=self.runPlot, width=bwidth
                   ).pack(side=tk.TOP)
    
    def buildMenus(self):
        '''
        builds the menu bar and contents
        '''
        if self.verbose: print("building the menus")
        # create a new menu
        self.menu = tk.Menu(self.root)

        # set the root menu to our new menu
        self.root.config(menu = self.menu)

        # create a variable to hold the top level menus
        menulist = []

        # create a file menu
        filemenu = tk.Menu( self.menu )
        self.menu.add_cascade( label = "File", menu = filemenu )
        menulist.append([filemenu,
                        [['Quit, Ctrl-Q OR Esc', self.selectImages], 
                         ['', None],
                         ['', None]
                         ]])

        # build the menu elements and callbacks
        for menu, items in menulist:
            for item in items:
                
                # blank menu item
                if not item[0]:
                    menu.add_separator()
                
                # sub cascade (TODO: could be recursive, only depth 1 now)
                elif isinstance(item[1], types.ListType):
                    submenu = tk.Menu(self.menu)
                    for subitem in item[1]:
                        submenu.add_command(label=subitem[0], command=subitem[1])
                    menu.add_cascade( label=item[0], menu=submenu )
                
                # menu command
                else:
                    menu.add_command( label=item[0], command=item[1] )
            
    def buildStatus(self):
        '''
        build the status frame
        '''
        self.statusFrame = tk.Frame(self.root)
        self.statusFrame.pack(side=tk.LEFT, padx=2, pady=2, fill=tk.Y)
        
        tk.Label(self.statusFrame, text="Status"
                 ).pack(side=tk.TOP)
        tk.Label(self.statusFrame, text="\nRun Directory:"
                 ).pack(side=tk.TOP)
        self.runDirectory = tk.StringVar()
        self.runDirectory.set("")
        tk.Label(self.statusFrame, textvariable=self.runDirectory
                 ).pack(side=tk.TOP)
                 
        tk.Label(self.statusFrame, text="\nImage Filenames:"
                 ).pack(side=tk.TOP)
        self.images = tk.StringVar()
        self.images.set("")
        tk.Label(self.statusFrame, textvariable=self.images
                 ).pack(side=tk.TOP)
                 
        tk.Label(self.statusFrame, text="\nResult Filenames:"
                 ).pack(side=tk.TOP)
        self.results = tk.StringVar()
        self.results.set("")
        tk.Label(self.statusFrame, textvariable=self.results
                 ).pack(side=tk.TOP)
                 
        tk.Label(self.statusFrame, text="\nSummary Filename:"
                 ).pack(side=tk.TOP)
        self.summary = tk.StringVar()
        self.summary.set("")
        tk.Label(self.statusFrame, textvariable=self.summary
                 ).pack(side=tk.TOP)
                 
    def clearStatus(self):
        '''
        clear the current status
        '''
        self.images.set("")
        self.runDirectory.set("")
        self.results.set("")
        self.summary.set("")
            
    def handleQuit(self, event=None):
        '''
        quit the display application
        '''
        if self.verbose: print('Terminating')
        self.root.destroy()
        
    def main(self):
        '''
        start the application
        '''
        if self.verbose: print('Entering main loop')
        self.root.mainloop()
        
    def runModel(self):
        '''
        run the modeling program
        '''
        if not self.images.get():
            tkm.showerror("No Filenames", "No filenames to model")
            return
        elif not self.runDirectory.get():
            tkm.showerror("No Directory", "No directory selected for modeling")
            return
        curWD = os.getcwd()
        modelPy = os.path.join(curWD, "simUtility.py")
        if not os.path.isfile(modelPy):
            tkm.showerror("Modeling program DNE", "Could not find %s" % modelPy)
            return
        if self.verbose: print("running the modeling program")
        os.chdir(self.runDirectory.get())
        imFilename = os.path.join(self.runDirectory.get(), "images.txt")
        with open(imFilename, "w") as imFile:
            imFile.write(self.images.get())
        commandList = ["python", modelPy, imFilename]
        if self.modelBulgeOpt.get(): commandList.append("-b")
        if self.modelGalfitOffOpt.get(): commandList.append("-g")
        if self.modelParallelOpt.get(): commandList.append("-p")
        if self.modelRealOpt.get(): commandList.append("-r")
        mpz = self.modelMPZEntry.get()
        try:
            float(mpz)
            commandList.append("--mpz")
            commandList.append(mpz)
        except ValueError:
            pass
        plate = self.modelPlateEntry.get()
        try:
            float(plate)
            commandList.append("--plate")
            commandList.append(plate)
        except ValueError:
            pass
        
        os.system(" ".join(commandList))
        if self.verbose: print("done modeling")
        os.chdir(curWD)
        
    def runSummary(self):
        '''
        run the summarizing program
        '''
        if self.verbose: print("running the summarizing program")
        if not self.results.get():
            tkm.showerror("No Filenames", "No results to model")
            return
        curWD = os.getcwd()
        sumPy = os.path.join(curWD, "sumSimUtility.py")
        if not os.path.isfile(sumPy):
            tkm.showerror("Summarizing program DNE", "Could not find %s" % sumPy)
            return
        resFilename = os.path.join(self.runDirectory.get(), "resultFilenames.txt")
        with open(resFilename, "w") as imFile:
            imFile.write(self.images.get())
        commandList = ["python", sumPy, resFilename]
        if self.verbose: commandList.append("-v")
        if self.sumBulgeOpt.get(): commandList.append("-b")
        if self.sumRealOpt.get(): commandList.append("-r")
        delim = self.sumDelimEntry.get()
        if delim: 
            commandList.append("-d")
            commandList.append(delim)
        output = self.sumOutputEntry.get()
        if delim: 
            commandList.append("-o")
            commandList.append(output)
        
        os.system(" ".join(commandList))
        if self.verbose: print("done summarizing")
        
    def runDisplay(self):
        '''
        run the displaying program
        '''
        if self.verbose: print("running the displaying program")
        
    def runPlot(self):
        '''
        run the plotting program
        '''
        if self.verbose: print("running the plotting program")

    def selectImages(self):
        '''
        select the fits images to be modeled
        '''
        if self.verbose: print("selecting images")
        filenames = tkf.askopenfilenames(parent=self.root, 
                                         title='Choose image files', 
                                         initialdir='.')
        if filenames:
            if self.verbose: print(filenames)
            self.images.set("\n".join(filenames))

    def selectResults(self):
        '''
        select the fits cube galfit results to be summarized
        '''
        if self.verbose: print("selecting results")
        filenames = tkf.askopenfilenames(parent=self.root, 
                                         title='Choose result files', 
                                         initialdir='.')
        if filenames:
            if self.verbose: print(filenames)
            self.results.set("\n".join(filenames))

    def selectRunModelDir(self):
        '''
        select the directory in which to run the modeling program
        '''
        if self.verbose: print("selecting run directory for modeling")
        directory = tkf.askdirectory(parent=self.root, 
                                     title='Choose Directory', 
                                     initialdir='.')
        if directory:
            if self.verbose: print(directory)
            self.runDirectory.set(directory)
        
    def selectSummary(self):
        '''
        select the csv summary file to display
        '''
        if self.verbose: print("selecting summary")
        filename = tkf.askopenfilename(parent=self.root, 
                                       title='Choose summary file', 
                                       initialdir='.')
        if filename:
            if self.verbose: print(filename)
            self.summary.set(filename)

    def setBindings(self):
        '''
        set the key bindings for the application
        '''
        if self.verbose: print("setting the key bindings")
        
        # bind command sequences to the root window
        self.root.bind( '<Control-q>', self.handleQuit)
        self.root.bind( '<Escape>', self.handleQuit)

if __name__ == "__main__":
    
    #define the command line interface usage message
    usage = "python %prog [-h help] [options (with '-'|'--' prefix)]"

    # used to parse command line arguments, -h will print options by default
    parser = OptionParser(usage)
    
    # indicate that the program should print actions to console
    parser.add_option("-v","--verbose", 
                help="to enable command line printouts of state",
                action="store_true")
    
    [options, args] = parser.parse_args()
    
    # run the application
    ModelerDashboard(1200, 800, options.verbose).main()