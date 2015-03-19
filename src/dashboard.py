#!/usr/bin/env python

'''
Author: Ian Tibbetts
Co-authors: Prof. Elizabeth McGrath
Last Edited: 8/13/2104
Colby College Astrophysics Research
'''

# use tkinter for GUI elements
import types
import threading
from optparse import OptionParser
import os
import simUtility
import sumSimUtility
# import plotSimUtility

try:
    import tkinter as tk  # python 3
    tkf = tk.filedialog
    tkm = tk.messagebox
    ttk = tk.ttk
except ImportError:
    import Tkinter as tk  # python 2
    import tkFileDialog as tkf
    import tkMessageBox as tkm
    import ttk
    
class Dashboard:
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
        self.root.geometry("%dx%d+50+30" % (width, height))
        
        # set the title of the window
        self.root.title("Modeling Dashboard")
        # set the maximum size of the window for resizing
        self.root.maxsize(1600, 900)
        
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
        headers = [header.capitalize() + "\n" for header in headers]
        frames = [tk.Frame(self.root) for _ in range(len(headers))]
        self.modelFrame, self.sumFrame, self.displayFrame, self.plotFrame = frames
        
        # pack all frames with separators and headers
        for header, frame in zip(headers, frames):
            frame.pack(side=tk.LEFT, padx=2, pady=2, fill=tk.Y)
            tk.Frame(self.root, width=2, bd=1, relief=tk.SUNKEN
                      ).pack(side=tk.LEFT, padx=2, pady=2, fill=tk.Y)
            tk.Label(frame, text=header
                      ).pack(side=tk.TOP)
                      
        bwidth = 20
        # make the model frame controls
        tk.Button(self.modelFrame, text="Select Images",
                   command=self.selectImages, width=bwidth
                   ).pack(side=tk.TOP)
        tk.Label(self.modelFrame, text="OR").pack(side=tk.TOP)
        tk.Button(self.modelFrame, text="Select Image File",
                   command=self.selectImageFile, width=bwidth
                   ).pack(side=tk.TOP)
        tk.Label(self.modelFrame, text="------ AND ------").pack(side=tk.TOP)
        tk.Button(self.modelFrame, text="Select Run Directory",
                   command=self.selectRunModelDir, width=bwidth
                   ).pack(side=tk.TOP)
        self.modelBulgeOpt = tk.IntVar()
        tk.Checkbutton(self.modelFrame, text="include bulge component",
                        variable=self.modelBulgeOpt
                        ).pack(side=tk.TOP)
        self.modelGalfitOffOpt = tk.IntVar()
        tk.Checkbutton(self.modelFrame, text="dont run GALFIT",
                        variable=self.modelGalfitOffOpt
                        ).pack(side=tk.TOP)
        self.modelParallelOpt = tk.IntVar()
        tk.Checkbutton(self.modelFrame, text="run in parallel",
                        variable=self.modelParallelOpt
                        ).pack(side=tk.TOP)
        self.modelRealOpt = tk.IntVar()
        tk.Checkbutton(self.modelFrame, text="images degraded",
                        variable=self.modelRealOpt
                        ).pack(side=tk.TOP)
        tk.Label(self.modelFrame, text="MPZ").pack(side=tk.TOP)
        self.modelMPZEntry = tk.Entry(self.modelFrame, width=bwidth)
        self.modelMPZEntry.pack(side=tk.TOP)
        tk.Label(self.modelFrame, text="Plate Scale").pack(side=tk.TOP)
        self.modelPlateEntry = tk.Entry(self.modelFrame, width=bwidth)
        self.modelPlateEntry.pack(side=tk.TOP)
        tk.Button(self.modelFrame, text="Run Modeling",
                   command=self.runModel, width=bwidth
                   ).pack(side=tk.TOP) 
        self.modelPBProgress = tk.IntVar()
        self.modelPB = ttk.Progressbar(master=self.modelFrame,
                                       orient="horizontal", 
                                       mode="determinate",
                                       maximum=0,
                                       variable=self.modelPBProgress)
        self.modelPB.pack(side=tk.TOP)
        #self.modelPBProgress.trace("w", self.updateRemaining)
        self.modelPBRemaining = tk.StringVar()
        tk.Label(self.modelFrame, textvariable=self.modelPBRemaining
                 ).pack(side=tk.TOP)
        self.modelPBProgress.set(0)
        
        # make the summary frame controls
        tk.Button(self.sumFrame, text="Select Results",
                   command=self.selectResults, width=bwidth
                   ).pack(side=tk.TOP)
        tk.Label(self.sumFrame, text="OR").pack(side=tk.TOP)
        tk.Button(self.sumFrame, text="Select Result File",
                   command=self.selectResultFile, width=bwidth
                   ).pack(side=tk.TOP)
        self.sumBulgeOpt = tk.IntVar()
        tk.Checkbutton(self.sumFrame, text="bulge components included",
                        variable=self.sumBulgeOpt
                        ).pack(side=tk.TOP)
        self.sumRealOpt = tk.IntVar()
        tk.Checkbutton(self.sumFrame, text="images degraded",
                        variable=self.sumRealOpt
                        ).pack(side=tk.TOP)
        tk.Label(self.sumFrame, text="Delimiter").pack(side=tk.TOP)
        self.sumDelimEntry = tk.Entry(self.sumFrame, width=bwidth)
        self.sumDelimEntry.pack(side=tk.TOP)
        tk.Label(self.sumFrame, text="Output Filename").pack(side=tk.TOP)
        self.sumOutputEntry = tk.Entry(self.sumFrame, width=bwidth)
        self.sumOutputEntry.pack(side=tk.TOP)
        tk.Button(self.sumFrame, text="Run Summarizing",
                   command=self.runSummary, width=bwidth
                   ).pack(side=tk.TOP)
        self.sumPBProgress = tk.IntVar()
        self.sumPB = ttk.Progressbar(master=self.sumFrame,
                                       orient="horizontal", 
                                       mode="determinate",
                                       maximum=0,
                                       variable=self.sumPBProgress)
        self.sumPB.pack(side=tk.TOP)
        #self.sumPBProgress.trace("w", self.updateRemaining)
        self.sumPBRemaining = tk.StringVar()
        tk.Label(self.sumFrame, textvariable=self.sumPBRemaining
                 ).pack(side=tk.TOP)
        self.sumPBProgress.set(0)
        
        # make the display frame controls
        tk.Button(self.displayFrame, text="Select Summary",
                   command=self.selectSummary, width=bwidth
                   ).pack(side=tk.TOP)
        tk.Button(self.displayFrame, text="Run Displaying",
                   command=self.runDisplay, width=bwidth
                   ).pack(side=tk.TOP)
        
        # make the plot frame controls
        tk.Button(self.plotFrame, text="Select Summary",
                   command=self.selectSummary, width=bwidth
                   ).pack(side=tk.TOP)
        tk.Button(self.plotFrame, text="Run Plotting",
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
        self.root.config(menu=self.menu)

        # create a variable to hold the top level menus
        menulist = []

        # create a file menu
        filemenu = tk.Menu(self.menu)
        self.menu.add_cascade(label="File", menu=filemenu)
        menulist.append([filemenu,
                        [['Quit, Ctrl-Q OR Esc', self.selectImages],
                         ['Reset', self.reset],
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
                    menu.add_cascade(label=item[0], menu=submenu)
                
                # menu command
                else:
                    menu.add_command(label=item[0], command=item[1])
            
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
                 
        tk.Label(self.statusFrame, text="\nImage Filename:"
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
        
    def reset(self):
        '''
        reset the application
        '''
        if self.verbose: print("resetting dashboard")
        self.images.set("")
        self.runDirectory.set("")
        self.results.set("")
        self.summary.set("")
        self.modelBulgeOpt.set(0)
        self.modelGalfitOffOpt.set(0)
        self.modelParallelOpt.set(0)
        self.modelRealOpt.set(0)
        self.modelMPZEntry.delete(0, tk.END)
        self.modelPlateEntry.delete(0, tk.END)
        self.sumBulgeOpt.set(0)
        self.sumRealOpt.set(0)
        self.sumDelimEntry.delete(0, tk.END)
        self.sumOutputEntry.delete(0, tk.END)
        
    def runDisplay(self):
        '''
        run the displaying program
        '''
        if self.verbose: print("running the displaying program")
        
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
        if self.verbose: print("running the modeling program")
        if self.images.get()[-5:] == ".fits":
            imFilename = os.path.join(self.runDirectory.get(), "images.txt")
            with open(imFilename, "w") as imFile:
                imFile.write(self.images.get())
        else:
            imFilename = self.images.get()
        commandList = [imFilename]
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
           
        if self.verbose: print(commandList)
        self.modelPB["maximum"] = self.numImages
        self.modelPBProgress.set(0)
        t1 = threading.Thread(target=simUtility.mainExternal, 
                              args=(commandList, self.runDirectory.get(), 
                                    os.getcwd(), self.modelPBProgress))
        t1.start()
        
#         tkm.showinfo("Started Modeling", 
# '''
# The modeling is in progress. 
# It will continue to run in a separate thread.
# The console will show its output, and a progress bar will show its progress.
# Each image takes around a minute, so it will take a while to finish running.
# You may quit the dashboard, but aborting the modeling may strand processes.
# ''')        
        
    def runPlot(self):
        '''
        run the plotting program
        '''
        if self.verbose: print("running the plotting program")
        self.modelPBProgress.set(self.modelPBProgress.get()+1)
        
    def runSummary(self):
        '''
        run the summarizing program
        '''
        if self.verbose: print("running the summarizing program")
        if not self.results.get():
            tkm.showerror("No Filenames", "No results to model")
            return
        if self.results.get()[-5:] == ".fits":
            resFilename = "resultFilenames.txt"
            with open(resFilename, "w") as imFile:
                imFile.write(self.results.get())
        else:
            resFilename = self.results.get()
        commandList = [resFilename]
        if self.verbose: commandList.append("-v")
        if self.sumBulgeOpt.get(): commandList.append("-b")
        if self.sumRealOpt.get(): commandList.append("-r")
        delim = self.sumDelimEntry.get()
        if delim: 
            commandList.append("-d")
            commandList.append(delim)
        output = self.sumOutputEntry.get()
        if output: 
            commandList.append("-o")
            commandList.append(output)
        if self.verbose: print(commandList)
        self.sumPB["maximum"] = self.numResults
        self.sumPBProgress.set(0)
        t1 = threading.Thread(target=sumSimUtility.main, 
                              args=(commandList, self.sumPBProgress))
        t1.start()

    def selectImages(self):
        '''
        select the fits images to be modeled
        '''
        if self.verbose: print("selecting images")
        if self.images.get():
            imageDir = self.images.get().split("\n")[0]
            imageDir = os.path.dirname(imageDir)
        else:
            imageDir = "."
        filenames = tkf.askopenfilenames(parent=self.root,
                                         title='Choose image files',
                                         initialdir=imageDir)
        if filenames:
            self.numImages = len(filenames)
            self.images.set("\n".join(filenames))
            
    def selectImageFile(self):
        '''
        select the file containing full path image files
        '''
        if self.verbose: print("selecting images file")
        if self.images.get() and self.images.get()[-5:] != ".fits":
            imageDir = os.path.dirname(self.images.get())
        else:
            imageDir = "."
        filename = tkf.askopenfilename(parent=self.root,
                                       title='Choose images file',
                                       initialdir=imageDir)
        if filename:
            with open(filename, 'r') as ifile:
                self.numImages = len(ifile.readlines())
            self.images.set(filename)

    def selectResults(self):
        '''
        select the fits cube galfit results to be summarized
        '''
        if self.verbose: print("selecting results")
        if self.results.get():
            resultDir = self.results.get().split("\n")[0]
            resultDir = os.path.dirname(resultDir)
        elif self.runDirectory.get():
            resultDir = self.runDirectory.get()
        else:
            resultDir = "."
        filenames = tkf.askopenfilenames(parent=self.root,
                                         title='Choose result files',
                                         initialdir=resultDir)
        if filenames:
            self.numResults = len(filenames)
            self.results.set("\n".join(filenames))
            
    def selectResultFile(self):
        '''
        select the file containing full path result files
        '''
        if self.verbose: print("selecting results file")
        if self.results.get() and self.results.get()[-5:] != ".fits":
            resultDir = os.path.dirname(self.results.get())
        else:
            resultDir = "."
        filename = tkf.askopenfilename(parent=self.root,
                                       title='Choose results file',
                                       initialdir=resultDir)
        if filename:
            with open(filename, 'r') as rfile:
                self.numResults = len(rfile.readlines())
            self.results.set(filename)

    def selectRunModelDir(self):
        '''
        select the directory in which to run the modeling program
        '''
        if self.runDirectory.get():
            curDir = self.runDirectory.get()
        else:
            curDir = "."
        if self.verbose: print("selecting run directory for modeling")
        directory = tkf.askdirectory(parent=self.root,
                                     title='Choose Directory',
                                     initialdir=curDir)
        if directory:
            if self.verbose: print(directory)
            self.runDirectory.set(directory)
        
    def selectSummary(self):
        '''
        select the csv summary file to display
        '''
        if self.verbose: print("selecting summary")
        if self.summary.get():
            sumDir = os.path.dirname(self.summary.get())
        else:
            sumDir = "."
        filename = tkf.askopenfilename(parent=self.root,
                                       title='Choose summary file',
                                       initialdir=sumDir)
        if filename:
            if self.verbose: print(filename)
            self.summary.set(filename)

    def setBindings(self):
        '''
        set the key bindings for the application
        '''
        if self.verbose: print("setting the key bindings")
        
        # bind command sequences to the root window
        self.root.bind('<Control-q>', self.handleQuit)
        self.root.bind('<Escape>', self.handleQuit)
        
    def updateRemaining(self, *args):
        '''
        update the estimated remaining time for the modeling
        '''
        m = self.modelPB["maximum"]-self.modelPBProgress.get()
        h, m = divmod(m, 60)
        self.modelPBRemaining.set("Remaining hr:min = "+str(h)+":"+str(m))

if __name__ == "__main__":
    
    # define the command line interface usage message
    usage = "python %prog [-h help] [options (with '-'|'--' prefix)]"

    # used to parse command line arguments, -h will print options by default
    parser = OptionParser(usage)
    
    # indicate that the program should print actions to console
    parser.add_option("-v", "--verbose",
                help="to enable command line printouts of state",
                action="store_true")
    
    [options, args] = parser.parse_args()
    
    # run the application
    md = Dashboard(1200, 500, options.verbose)
    md.main()
    print("done")
