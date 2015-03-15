import wx, os, glob
import pdb
import numpy as np
import matplotlib
import processCCD_image

# comment out the following to use wx rather than wxagg
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

# Add display box showing header of file?

#correct zoom using rubbish/rectange_select.py code
#self.axes.remove_patch(self.rect)
#somehow remove patch?

class ControlFrame(wx.Frame):
    def __init__(self):
        wx.Frame.__init__(self, None, title='File Browser', size=(600, 500), pos=(50, 50))

        panel = wx.Panel(self, -1)
        
        vbox = wx.BoxSizer(wx.VERTICAL)
        bdr_sz = 8 # was 20
        
        # path TextCtrl
        vbox.AddSpacer(bdr_sz)
        path_cap = wx.StaticText( panel, label='Search Term' )
        vbox.Add(path_cap, 0, wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, border=bdr_sz)
        self.path = wx.TextCtrl(panel, size=(140, 20), value = 
                      '/Users/mdean/Dropbox/SLS_March_2015/RIXS/*.txt',
                      style=wx.TE_PROCESS_ENTER)
        self.path.Bind(wx.EVT_TEXT_ENTER, self.updatepath)
        vbox.Add(self.path, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, border=bdr_sz)
        vbox.AddSpacer(bdr_sz)

#        # eV per pix TextCtrl
#        #vbox.AddSpacer(bdr_sz)
#        eVperPix_cap = wx.StaticText( panel, label='eV per Pixel' )
#        vbox.Add(eVperPix_cap, 0, wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, border=bdr_sz)
#        self.eVperPix = wx.TextCtrl(panel, size=(140, 20), value = 
#                      '0.0771', style=wx.TE_PROCESS_ENTER)
#        #self.path.Bind(wx.EVT_TEXT_ENTER, self.updateeVperPix)
#        vbox.Add(self.eVperPix, 0, wx.LEFT | wx.RIGHT | wx.EXPAND, border=bdr_sz)
#        vbox.AddSpacer(bdr_sz)

        
        # LIST BOX SHOWING FILES
        self.path_list, self.file_list = self.get_files(self.path.GetValue())

        listbox_cap = wx.StaticText( panel, label='Directory contents' )
        vbox.Add(listbox_cap, 0, wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, border=bdr_sz)
        self.listbox = wx.ListBox(panel, -1, wx.DefaultPosition, (170, 130), self.file_list, wx.LB_MULTIPLE)
        
        
        # buttons
        vbox.Add(self.listbox, 1, wx.EXPAND | wx.LEFT |wx.RIGHT | wx.BOTTOM, border=bdr_sz)
        buttonsize=(80, 20)
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        self.plot_but = wx.Button(panel, -1, 'Plot', size=buttonsize)
        self.Bind(wx.EVT_BUTTON, self.plotit, self.plot_but)
        hbox.Add(self.plot_but, 0, wx.ALIGN_LEFT, border=bdr_sz)
        hbox.AddSpacer(bdr_sz)
        
        self.clear_but = wx.Button(panel, -1, 'Clear', size=buttonsize)
        self.Bind(wx.EVT_BUTTON, self.clearit, self.clear_but)
        hbox.Add(self.clear_but, 0, wx.ALIGN_LEFT, border=bdr_sz)
        hbox.AddSpacer(bdr_sz)
        
        self.autoscale_but = wx.Button(panel, -1, 'Autoscale', size=buttonsize)
        self.Bind(wx.EVT_BUTTON, self.autoscaleit, self.autoscale_but)
        hbox.Add(self.autoscale_but, 0, wx.ALIGN_LEFT, border=bdr_sz)
        hbox.AddSpacer(bdr_sz)
        
        ##### NEW
        #hbox = wx.BoxSizer(wx.HORIZONTAL)

        self.align_but = wx.Button(panel, -1, 'Align', size=buttonsize)
        self.Bind(wx.EVT_BUTTON, self.alignit, self.align_but)
        hbox.Add(self.align_but, 0, wx.ALIGN_LEFT, border=bdr_sz)
        hbox.AddSpacer(bdr_sz)       

        self.sum_but = wx.Button(panel, -1, 'Sum', size=buttonsize)
        self.Bind(wx.EVT_BUTTON, self.sumit, self.sum_but)
        hbox.Add(self.sum_but, 0, wx.ALIGN_LEFT, border=bdr_sz)
        hbox.AddSpacer(bdr_sz)  
        
        self.calibrate_but = wx.Button(panel, -1, 'Calibrate', size=buttonsize)
        self.Bind(wx.EVT_BUTTON, self.calibrateit, self.calibrate_but)
        hbox.Add(self.calibrate_but, 0, wx.ALIGN_LEFT, border=bdr_sz)
        hbox.AddSpacer(bdr_sz)  

        self.save_but = wx.Button(panel, -1, 'Save', size=buttonsize)
        self.Bind(wx.EVT_BUTTON, self.saveit, self.save_but)
        hbox.Add(self.save_but, 0, wx.ALIGN_LEFT, border=bdr_sz)
        hbox.AddSpacer(bdr_sz) 

#        self.legend_but = wx.Button(panel, -1, 'Legend', size=buttonsize)
#        self.Bind(wx.EVT_BUTTON, self.legendtoggleit, self.legend_but)
#        hbox.Add(self.legend_but, 0, wx.ALIGN_LEFT, border=bdr_sz)
        
        vbox.Add(hbox, 0, wx.LEFT | wx.RIGHT | wx.ALIGN_LEFT, border=bdr_sz)
        vbox.AddSpacer(bdr_sz)
        
        self.legend_exists = False
        
        # Initiate everything
        panel.SetSizer(vbox)
        #self.Centre()
        self.Show(True)

        # Call plotframe        
        self.plotframe = PlotFrame(self)
        self.plotframe.Show()

    def updatepath(self, evt):
        print "update path"
        self.path_list, self.file_list = self.get_files(self.path.GetValue())
        self.listbox.Clear()
        #self.listbox.InsertItems(self.file_list, 0)
        for filelist in self.file_list:
            self.listbox.Append(filelist)
        #pdb.set_trace()
    
    def plotit(self, evt):
        print "plotting"
        self.plotframe.figure.axes[0].cla()
        self.ims = processCCD_image.CCD()
        self.ims.specs = []
        self.ims.BGspecs = []
        self.ims.filenames = []
        self.ims.curvature = 'Set in SAXES control program'
        for index in self.listbox.GetSelections():
            M = np.loadtxt(self.path_list[index], skiprows=1)
            filename = os.path.split(self.path_list[index])[1]
            self.ims.specs.append([M[:,0], M[:,2]])
            self.ims.BGspecs.append([M[:,0]* 0.0, M[:,2]*0.0])
            self.ims.filenames.append(filename)
            self.plotframe.figure.axes[0].plot(M[:,0],M[:,2], label=filename)
        self.legendobj = self.plotframe.axes.legend(prop={'size':9})     
        self.plotframe.figure.canvas.draw()
        
    def clearit(self, evt):
        print "Clear button pressed"
        self.ims.specs = []
        self.plotframe.figure.axes[0].cla()
        self.plotframe.figure.canvas.draw()

    def autoscaleit(self, evt):
        print "Autoscale pressed"
        self.plotframe.figure.axes[0].autoscale(True)
        self.plotframe.figure.canvas.draw()

    def alignit(self, evt):
        print "Align pressed"
        self.ims.correlate_specs()
        self.plotframe.figure.axes[0].cla()
        for M_list, filename in zip(self.ims.specs, self.ims.filenames):
            self.plotframe.figure.axes[0].plot(M_list[0], M_list[1], label=filename)    
        self.legendobj = self.plotframe.axes.legend(prop={'size':9})
        self.plotframe.figure.canvas.draw()

    def sumit(self, evt):
        print "Sum pressed"
        self.ims.sum_specs()
        self.plotframe.figure.axes[0].plot(self.ims.spectrum[0], self.ims.spectrum[1], label='Sum')
        self.legendobj = self.plotframe.axes.legend(prop={'size':9}) 
        self.plotframe.figure.canvas.draw()

    def calibrateit(self, evt):
        print "Calibrate pressed"
        dlg = wx.TextEntryDialog(None, "Enter eV per pixel", 'Calibration', "0.0771", 
                style=wx.OK)
        dlg.ShowModal()
        eVperPixel = float(dlg.GetValue())
        dlg = wx.TextEntryDialog(None, "Enter elastic pixel", 'Calibration', "", 
                style=wx.OK)
        dlg.ShowModal()
        elPixel = float(dlg.GetValue())
        dlg.Destroy()
        self.ims.calibrate(elPixel, eVperPixel)
        self.plotframe.figure.axes[0].cla()
        self.plotframe.figure.axes[0].plot(self.ims.spectrum[0], self.ims.spectrum[1], label='Calibrated sum')
        self.legendobj = self.plotframe.axes.legend(prop={'size':9}) 
        self.plotframe.figure.canvas.draw()
        
        
#        self.ims.sum_specs()
#        self.plotframe.figure.axes[0].plot(self.ims.spectrum[0], self.ims.spectrum[1], label='Sum')
#        self.legendobj = self.plotframe.axes.legend(prop={'size':9}) 
#        self.plotframe.figure.canvas.draw()
    
    def saveit(self, evt):
        print "Save pressed"
        saveFileDialog = wx.FileDialog(self, "Save file", "", "",
                                   "TXT files (*.txt)|*.txt", wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

        if saveFileDialog.ShowModal() == wx.ID_CANCEL:
            return     # the user changed idea...

        print saveFileDialog.GetPath()
        self.ims.fileout = saveFileDialog.GetPath()
        self.ims.write_file()


    def legendtoggleit(self, evt):
        print "Legend pressed"
        self.legendobj = self.plotframe.axes.legend(prop={'size':9})
        if self.legend_exists:
            self.legendobj.set_visible(False)
            self.legend_exists = False
        else:
            self.legendobj.set_visible(True)
            self.legend_exists = True

        self.plotframe.figure.canvas.draw()

#        self.legendobj = self.plotframe.axes.legend()
#        try:
#            self.legendobj.set_visible(False)
#        except AttributeError:
#            self.legendobj = self.plotframe.axes.legend()
        

    def onselect(vmin, vmax):
        print "onselect ran"
        print vmin, vmax

    def get_files(self, searchterm):
        path_list = glob.glob(searchterm)
        
        file_list = []
        for path in path_list:
            file_list.append(os.path.split(path)[1])
            
        return path_list, file_list
    
        

class PlotFrame(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(self, None, size=(300,500), pos=(650, 50), title='Plot Frame')
        self.parent = parent
        
        # initialize plot        
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)        
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        
        # Connect the mouse events to their relevant callbacks
        self.canvas.mpl_connect('button_press_event', self._onPress)
        self.canvas.mpl_connect('button_release_event', self._onRelease)
        self.canvas.mpl_connect('motion_notify_event', self._onMotion)
        
        self.pressed = False
        
        # Initialise the rectangle
        self.rect = Rectangle((0,0), 1, 1, facecolor='None', visible=False,
                                  edgecolor='k', linestyle='dashed')
        self.x0 = 0#None
        self.y0 = 0#None
        self.x1 = 0#None
        self.y1 = 0#None
        self.axes.add_patch(self.rect)
        
        self.Fit()
        
    def _onPress(self, event):
        ''' Callback to handle the mouse being clicked and held over the canvas'''
        print "onPress"
        # Check the mouse press was actually on the canvas
        if event.xdata is not None and event.ydata is not None:
            print "In press event data"
            # Upon initial press of the mouse record the origin and record the mouse as pressed
            self.pressed = True
            self.rect.set_visible = True
            self.x0 = event.xdata
            self.y0 = event.ydata
            
            # Set the width and height and draw the rectangle
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            #self.axes.add_patch(self.rect)
            self.canvas.draw()


    def _onRelease(self, event):
        '''Callback to handle the mouse being released over the canvas'''
        print "onRelease"
        # Check that the mouse was actually pressed on the canvas to begin with and this isn't a rouge mouse 
        # release event that started somewhere else
        if self.pressed:

            # Upon release draw the rectangle as a solid rectangle
            self.pressed = False

            # Check the mouse was released on the canvas, and if it wasn't then just leave the width and 
            # height as the last values set by the motion event
            if event.xdata is not None and event.ydata is not None:
                self.x1 = event.xdata
                self.y1 = event.ydata

            # Set the width and height and origin of the bounding rectangle
            self.boundingRectWidth =  self.x1 - self.x0
            self.boundingRectHeight =  self.y1 - self.y0
            self.bouningRectOrigin = (self.x0, self.y0)

            # Draw the bounding rectangle
            self.rect.set_width(self.boundingRectWidth)
            self.rect.set_height(self.boundingRectHeight)
            self.rect.set_xy((self.x0, self.y0))
            
            # set zoom
            xmin = min([self.x0, self.x1])
            xmax = max([self.x0, self.x1])
            ymin = min([self.y0, self.y1])
            ymax = max([self.y0, self.y1])
            self.figure.axes[0].axis([xmin, xmax, ymin, ymax])
            self.rect.set_visible = False

            self.canvas.draw()


    def _onMotion(self, event):
        '''Callback to handle the motion event created by the mouse moving over the canvas'''
        #print "onMotion"
        # If the mouse has been pressed draw an updated rectangle when the mouse is moved so 
        # the user can see what the current selection is
        if self.pressed:

            # Check the mouse was released on the canvas, and if it wasn't then just leave the width and 
            # height as the last values set by the motion event
            if event.xdata is not None and event.ydata is not None:
                self.x1 = event.xdata
                self.y1 = event.ydata
            
            # Set the width and height and draw the rectangle
            self.rect.set_visible = True
            self.rect.set_width(self.x1 - self.x0)
            self.rect.set_height(self.y1 - self.y0)
            self.rect.set_xy((self.x0, self.y0))
            #self.rectpatch = self.axes.add_patch(self.rect)
            self.canvas.draw()
        
    def OnPaint(self, event):
        self.canvas.draw()

if __name__ == "__main__":
    App=wx.PySimpleApp()
    ControlFrame().Show()
    App.MainLoop()