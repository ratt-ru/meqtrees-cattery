#!/usr/bin/python

if __name__ == "__main__":
  print "Welcome to SkyPuss!"
  print "Please wait a second while the GUI starts up."

  # parse options is the first thing we should do
  from optparse import OptionParser
  usage = "usage: %prog [options] <directories to watch>"
  parser = OptionParser()
  parser.add_option("-d", "--debug",dest="verbose",type="string",action="append",metavar="Context=Level",
                    help="(for debugging Python code) sets verbosity level of the named Python context. May be used multiple times.");
  (options, rem_args) = parser.parse_args();
  
  import sys
  import signal
  import os.path
  
  from qt import *
  
  import SkyPuss
  import SkyPuss.MainWindow
  import Timba.utils
  
  
  app = QApplication(sys.argv);
  app.setDesktopSettingsAware(True);
  
  mainwin = SkyPuss.MainWindow.MainWindow(None);
  app.setMainWidget(mainwin);
  mainwin.show();
  QObject.connect(app,SIGNAL("lastWindowClosed()"),app,SLOT("quit()"));
  
  # load initial file, if specified
  if rem_args:
    mainwin.openFile(rem_args[0]);
    
  # handle SIGINT
  def sigint_handler (sig,stackframe):
      mainwin.close();
  signal.signal(signal.SIGINT,sigint_handler);
      
  app.exec_loop(); 
  print "SkyPuss exiting normally, goodbye!";
