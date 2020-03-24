import sys
import time
import tkinter as tk
import pandas as pd
import networkx as nx
import neurontree.NeuronTree as nt
import matplotlib.pyplot as plt
from tkinter import messagebox, filedialog, ttk
import threading
import computation.persistence_functions as pf
import computation.feature_presentation as fr

class StdoutRedirector(object):
    def __init__(self, text_widget):
        self.text_space = text_widget

    def write(self, string):
        self.text_space.insert('end', string)
        self.text_space.see('end')

    def flush(self):
        #do nothing
        time.sleep(0)


# class for the main window
class MainWindow(object):
    def __init__(self, parent):
        self.parent = parent
        parent.title("Neuron Tree")
        parent.iconbitmap('@gui/neuron_tree.xbm')
        # object to access NeuronTree class
        self.myNeuronTree = nt.NeuronTree()
        self.figure = plt.figure()
        self.figure.canvas.mpl_connect('close_event', self.handle_close)
        # state if file is already loaded or not for showing more menu options
        self.fileloaded = False

        ## Center window on screen:
        # Gets both half the screen width/height and window width/height
        positionRight = int(parent.winfo_screenwidth() / 2 - 600 / 2)
        positionDown = int(parent.winfo_screenheight() / 2 - 480 / 2)
        # Positions the window in the center of the page.
        parent.geometry("+{}+{}".format(positionRight, positionDown))

        # frames for organizing GUI items
        frame = tk.Frame(parent)
        frame.pack(fill=tk.BOTH, expand=True)
        topframe = tk.Frame(frame)
        topframe.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        bottomframe = tk.Frame(frame)
        bottomframe.pack(side=tk.TOP, fill=tk.X)

        # add all GUI elements
        self.textbox = tk.Text(topframe)
        self.scrollbar = tk.Scrollbar(topframe)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.textbox.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.scrollbar.config(command=self.textbox.yview)
        self.textbox.config(yscrollcommand=self.scrollbar.set)
        sys.stdout = StdoutRedirector(self.textbox)
        sys.stderr = StdoutRedirector(self.textbox)

        self.progressbar = ttk.Progressbar(bottomframe, orient="horizontal", mode="indeterminate")
        self.progressbar.pack(fill=tk.X)

        self.menu = tk.Menu(frame)
        self.parent.config(menu=self.menu)

        self.filemenu = tk.Menu(self.menu)
        self.menu.add_cascade(label="File", menu=self.filemenu)
        self.filemenu.add_command(label="Open...", command=self.openfile)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Save to swc...", command=self.export2swc)
        self.filemenu.add_command(label="Save to mat...", command=self.export2mat)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=self.quitbutton)
        # X from Main Window redirect to same method
        self.parent.protocol("WM_DELETE_WINDOW", self.quitbutton)

        self.analysemenu = tk.Menu(self.menu)
        self.menu.add_cascade(label="Analyse", menu=self.analysemenu)
        self.analysemenu.add_command(label="Get Minimal Spanning Tree", command=self.minSpan)
        self.analysemenu.add_command(label="Compute Persistence...", command=self.createPersistence)
        self.analysemenu.add_command(label="Compute Density Map", command=self.computeDensity)
        self.analysemenu.add_command(label="Compute Morphpometric Statistics", command=self.computeMorphoStats)

        self.viewmenu = tk.Menu(self.menu)
        self.menu.add_cascade(label="View", menu=self.viewmenu)
        self.viewmenu.add_command(label="Draw tree ...", command=self.drawtree)
        self.viewmenu.add_command(label="Draw 2D ...", command=self.draw2d)
        self.viewmenu.add_command(label="Draw 3D ...", command=self.draw3d)
        self.viewmenu.add_command(label="Draw 3D volumetric ...", command=self.draw3dvol)

        # until files are not loaded disable menu
        self.filemenu.entryconfig("Save to swc...", state="disabled")
        self.filemenu.entryconfig("Save to mat...", state="disabled")
        self.menu.entryconfig("Analyse", state="disabled")
        self.menu.entryconfig("View", state="disabled")

        self.helpmenu = tk.Menu(self.menu)
        self.menu.add_cascade(label="Help", menu=self.helpmenu)
        self.helpmenu.add_command(label="About...", command=self.about)

        sys.stdout = StdoutRedirector(self.textbox)

        # check for required libraries
        print("Pandas availabe: Version " + pd.__version__)
        print("NetworkX available: Version "+ nx.__version__)
        try:
            import pygraphviz as pgv
            print("Pygraphviz availabe: Version " + pgv.__version__)
        except:
            print("Pygraphviz not installed.")

        # set version of networkX
        if float(nx.__version__) < 2:
            self._nxversion = 1
        else:
            self._nxversion = 2

    def enable_menu(self):
        self.menu.entryconfig("File", state="normal")
        if self.fileloaded:
            self.menu.entryconfig("Analyse", state="normal")
            self.menu.entryconfig("View", state="normal")
            self.filemenu.entryconfig("Save to swc...", state="normal")
            self.filemenu.entryconfig("Save to mat...", state="normal")

    def disable_menu(self):
        self.menu.entryconfig("File", state="disabled")
        self.menu.entryconfig("Analyse", state="disabled")
        self.menu.entryconfig("View", state="disabled")

    def openfile(self):
        print("Open file ...")
        self.filename = filedialog.askopenfilename(initialdir=".", title="Select file",
                                                   filetypes=(("swc files", "*.swc"), ("all files", "*.*")))
        if len(self.filename) > 0:
            self.disable_menu()
            swc = pd.read_csv(self.filename, delim_whitespace=True, comment='#',
                          names=['n', 'type', 'x', 'y', 'z', 'radius', 'parent'], index_col=False)
            thread = threading.Thread(target=self.loadThread, args=(swc, ))
            thread.start()
            self.progressbar.start()

    def loadThread(self, swc):
        try:
            self.myNeuronTree = nt.NeuronTree(swc=swc, nxversion=self._nxversion)
            self.fileloaded = True
            print(self.filename + " successfully loaded.")
        except:
            raise
            print("Fehler beim Laden!")
            print(sys.exc_info()[0])
        self.progressbar.stop()
        self.enable_menu()

    def handle_close(evt):
        plt.close('all')
        print('Closed Figure!')

    def quitbutton(self):
        if messagebox.askyesno('Verify', 'Really quit?'):
            plt.close('all')
            self.parent.destroy()

    def minSpan(self):
        print("Get minimal spanning tree ...")
        min_span_tree = self.myNeuronTree.get_topological_minor()

    def createPersistence(self):

        result = PersistenceDialog(self.parent).show()
        if result > 0:
            self.function = getattr(pf, pf.functions[result])
            print("Computes the persistence barcode with {} ...".format(pf.functions[result]))
        else:
            self.function = None
            print("Computes the persistence barcode without user function ...")

        self.figure = plt.figure()
        self.thread = threading.Thread(target=self.createPersThread)
        self.thread.start()
        self.disable_menu()
        self.parent.after(100, self.onDrawPlot)
        self.progressbar.start()

    def createPersThread(self):
        try:
            min_span_tree = self.myNeuronTree.get_topological_minor()
            persistence_table = fr.get_persistence(min_span_tree, self.function)
            print(persistence_table)
            x = persistence_table['birth']
            y = persistence_table['death']
            plt.scatter(x, y, alpha=0.5)
            if self.function is None:
                plt.title('Persistence Diagram')
            else:
                plt.title('Persistence Diagram ({})'.format(self.function.__name__))
            plt.xlabel('birth')
            plt.ylabel('death')
        except:
            print("Failure during computing persistence!")
            print(sys.exc_info()[0])

    def computeDensity(self):
        print("Computing a density map ...")
        self.figure = plt.figure()
        self.thread = threading.Thread(target=self.computeDensityThread)
        self.thread.start()
        self.disable_menu()
        self.parent.after(100, self.onDrawPlot)
        self.progressbar.start()

    def computeDensityThread(self):
        try:
            dist = 1 # in microns
            # get the resampled point could along each neurite at distance 1 micron.
            # pc is an array of 3D coordinates for each resampled node
            pc = self.myNeuronTree.resample_nodes(self.myNeuronTree.get_graph(), dist)
            plt.scatter(pc[:,0], pc[:,2], s=1)
            plt.title('Density Map')
        except:
            print("Failure during computing density map!")
            print(sys.exc_info()[0])

    def computeMorphoStats(self):
        print("Computing morphometric statistics ...")
        self.thread = threading.Thread(target=self.computeMorphoStatsThread)
        self.thread.start()
        self.disable_menu()
        self.progressbar.start()

    def computeMorphoStatsThread(self):
        try:
            results = fr.compute_Morphometric_Statistics(self.myNeuronTree)
            print(results)
        except:
            print("Failure during computing statistics!")
            print(sys.exc_info()[0])
        self.progressbar.stop()
        self.enable_menu()


    def onDrawPlot(self):
        if (self.thread.is_alive()):
            self.parent.after(100, self.onDrawPlot)
            return
        else:
            self.progressbar.stop()
            self.enable_menu()
            self.figure.show()

    def drawtree(self):
        print("Draw tree ...")
        self.disable_menu()
        self.progressbar.start()
        thread = threading.Thread(target=self.drawTreeThread)
        thread.start()

    def drawTreeThread(self):
        try:
            self.myNeuronTree.draw_tree()
            print("Tree successfully generated.")
        except:
            print("Failure during drawing of the tree!")
            print(sys.exc_info()[0])

    def draw2d(self):
        print("Draw 2D ...")
        self.figure = plt.figure()
        self.thread = threading.Thread(target=self.draw2dThread)
        self.thread.start()
        self.disable_menu()
        self.parent.after(100, self.onDrawPlot)
        self.progressbar.start()

    def draw2dThread(self):
        try:
            self.myNeuronTree.draw_2D(fig=self.figure)
            print("2D tree successfully generated.")
        except:
            print("Failure during drawing 2D!")
            print(sys.exc_info()[0])

    def draw3d(self):
        print("Draw 3D ...")
        self.figure = plt.figure()
        self.thread = threading.Thread(target=self.draw3dThread)
        self.thread.start()
        self.disable_menu()
        self.parent.after(100, self.onDrawPlot)
        self.progressbar.start()

    def draw3dThread(self):
        try:
            self.myNeuronTree.draw_3D(fig=self.figure)
            print("3D tree successfully generated.")
        except:
            print("Failure during drawing 3D!")
            print(sys.exc_info()[0])

    def draw3dvol(self):
        print("Draw 3D volumetric...")
        self.figure = plt.figure()
        self.thread = threading.Thread(target=self.draw3dvolThread)
        self.thread.start()
        self.disable_menu()
        self.parent.after(100, self.onDrawPlot)
        self.progressbar.start()

    def draw3dvolThread(self):
        try:
            self.myNeuronTree.draw_3D_volumetric(fig=self.figure)
            print("3D vol tree successfully generated.")
        except:
            print("Failure during drawing of 3D volumetric!")
            print(sys.exc_info()[0])
            raise

    def export2swc(self):
        print("Export to swc ...")
        self.filename = filedialog.asksaveasfilename(initialdir=".", title="Select file",
                                                     filetypes=(("swc files","*.swc"), ("all files", "*.*")))
        idx = self.filename.rfind("/")+1
        fname = self.filename[idx:]
        if fname.endswith('.swc'):
            fname = fname[:-4]
        path = self.filename[:idx]
        self.myNeuronTree.write_to_swc(fname, path)
        print("Save to (swc file): " + self.filename)

    def export2mat(self):
        print("Export to mat ...")
        self.filename = filedialog.asksaveasfilename(initialdir=".", title="Select file",
                                                     filetypes=(("mat files", "*.mat"), ("all files", "*.*")))
        idx = self.filename.rfind("/")+1
        fname = self.filename[idx:]
        if fname.endswith('.mat'):
            fname = fname[:-4]
        path = self.filename[:idx]
        self.myNeuronTree.write_to_mat(fname, path)
        print("Save to (matlab file): " + self.filename)

    def about(self):
        messagebox.showinfo("About", "        (c) 2019\nAdam von Daranyi")


class PersistenceDialog(tk.Toplevel):
    def __init__(self, parent):
        tk.Toplevel.__init__(self, parent)
        self.v = tk.IntVar()
        self.v.set(0)
        self.grab_set()
        self.overrideredirect(1)

        ## Center window on screen:
        # Gets both half the screen width/height and window width/height
        positionRight = int(parent.winfo_screenwidth() / 2 - 100 / 2)
        positionDown = int(parent.winfo_screenheight() / 2 - 200 / 2)
        # Positions the window in the center of the page.
        self.geometry("+{}+{}".format(positionRight, positionDown))

        tk.Label(self, text="""Choose a function:""", justify=tk.LEFT, padx=20).pack()
        for i, func in enumerate(pf.functions):
            tk.Radiobutton(self, text=func, padx=20, variable=self.v, val=i).pack(anchor=tk.W)

        self.ok_button = tk.Button(self, text="OK", command=self.on_ok).pack(side="right")

    def on_ok(self, event=None):
        self.destroy()

    def show(self):
        self.wm_deiconify()
        self.wait_window()
        return self.v.get()
