from ModelClasses import *
import PlotStyles
 
class ModelTag (ModelItem):
  mandatory_attrs = [ "name" ];
  optional_attrs = dict([ (attr,None) for attr in PlotStyles.StyleAttributes ]);
 
ModelTag.registerClass();
  
class ModelTagSet (ModelItem):
  def __init__ (self,*tags,**kws):
    ModelItem.__init__(self,**kws);
    self.tags = dict([ (tag.name,tag) for tag in tags ]);
    
  def add (self,tag):
    """Adds a ModelTag object to the tag set""";
    self.tags[tag.name] = tag;
    
  def get (self,tagname):
    """Returns ModelTag object associated with tag name, inserting a new one if not found""";
    return self.tags.setdefault(name,ModelTag(name));
    
  def getAll (self):
    all = self.tags.values();
    all.sort(lambda a,b:cmp(a.name,b.name));
    return all;
    
  def addNames (self,names):
    """Ensures that ModelTag objects are initialized for all tagnames in names""";
    for name in names:
      self.tags.setdefault(name,ModelTag(name));
    
  def renderMarkup (self,tag="A",attrname=None):
      """Makes a markup string corresponding to the model item.
      'tags' is the HTML tag to use.
      """;
      # opening tag
      markup = "<%s mdltype=ModelTagList "%tag;
      if attrname is not None:
        markup += "mdlattr=%s "%attrname;
      markup +=">";
      # write mandatory attributes
      for name,tt in self.tags.iteritems():
        markup += self.renderAttrMarkup(name,tt,tag="TR",mandatory=True);
      # closing tag
      markup += "</%s>"%tag;
      return markup;
ModelTagSet.registerClass();
 
class Source (ModelItem):
  mandatory_attrs  = [ "name","pos","flux" ];
  optional_attrs     = dict(spectrum=None,shape=None);
  allow_extra_attrs = True;
  
  def __init__ (self,*args,**kw):
    ModelItem.__init__(self,*args,**kw);
    self.typecode = (self.shape and self.shape.typecode) or "pnt";
    self.selected = False;
    
  def select (self,sel=True):
    self.selected = sel;
    
  def brightness (self):
    iapp = getattr(self,'Iapp',None);
    if iapp is not None:
      return iapp;
    else:
      return getattr(self.flux,'I',0.);
      
  def getTagNames (self):
    return [ attr for attr,val in self.getExtraAttributes() if attr[0] != "_" ];
    
  def getTags (self):
    return [ (attr,val) for attr,val in self.getExtraAttributes() if attr[0] != "_" ];
    
  class Property (object):
    def __init__ (self,name,func,style=PlotStyles.DefaultPlotStyle,sources=None):
      self.name = name;
      self.style = style.copy();
      self.func = func;
      self.total = 0;
      self.show_plot = self.show_list = False;
      self.apply = True;
      if sources:
        self.computeTotal(sources);
    def computeTotal (self,sources):
      self.total = len(filter(self.func,sources));
      return self.total;

Source.registerClass();

class SkyModel (ModelItem):
  optional_attrs   = dict(name=None,plotstyles={},pbexp=None,ra0=None,dec0=None);
  
  def __init__ (self,*sources,**kws):
    ModelItem.__init__(self,**kws);
    self.updated = False;
    self.sources = list(sources);
    self.plotstyles = self.plotstyles or {};
    # build up set of source tags
    self.scanTags();
    self.initProperties();
  
  def add (self,*src):
    self.sources += list(src);
    self._scanTags(src);
    self._computeProperties();
    
  def scanTags (self,sources=None):
    """Populates self.tagnames with a list of tags present in sources""";
    sources = sources or self.sources;
    tagnames = sets.Set();
    for src in sources:
      tagnames.update(src.getTagNames());
    self.tagnames = list(tagnames);
    self.tagnames.sort();
    
  def initProperties (self):
    # init default and "selected" properties
    self.defprop = Source.Property("all sources",func=lambda src:True,sources=self.sources,
                                                        style=self.plotstyles.setdefault('default',PlotStyles.BaselinePlotStyle));
    self.defprop.show_plot = self.defprop.show_list = True;
    self.selprop = Source.Property("selected",func=lambda src:getattr(src,'selected',False),sources=self.sources,
                                                        style=self.plotstyles.setdefault('selected',PlotStyles.SelectionPlotStyle));
    self.selprop.show_plot = True; 
    self.selprop.show_list = True,"always";  # ModelTagsTable will make this uneditable
    # and make ordered list of properties
    self.properties = [ self.defprop,self.selprop ];
    # make properties from available source types
    self._typeprops = {};
    typecodes = list(sets.Set([src.typecode for src in self.sources]));
    typecodes.sort();
    if len(typecodes) > 1:
      for code in typecodes:
          self._typeprops[code] = prop = Source.Property("type: %s"%code,lambda src,code=code:src.typecode==code,sources=self.sources,
                                                                            style=self.plotstyles.setdefault('type:%s'%code,PlotStyles.DefaultPlotStyle));
          self.properties.append(prop);
    # make properties from source tags
    self._tagprops = {};
    for tag in self.tagnames:
      self._tagprops[tag] = prop = Source.Property("tag: %s"%tag,
                                                                    lambda src,tag=tag:getattr(src,tag,None) not in [None,False],
                                                                    sources=self.sources,
                                                                    style=self.plotstyles.setdefault('tag:%s'%tag,PlotStyles.DefaultPlotStyle));
      self.properties.append(prop);
      
  def _remakePropList (self):
    self.properties = [ self.defprop,self.selprop ];
    typenames = self._typeprops.keys();
    typenames.sort();
    self.properties += [ self._typeprops[name] for name in typenames ];
    self.properties += [ self._tagprops[name] for name in self.tagnames ];
    
  def getTagProperty (self,tag):
    return self._tagprops[tag];
    
  def getTypeProperty (self,typename):
    return self._typeprops[typename];
  
  def addTag (self,tag):
    if tag in self.tagnames:
      return False;
    # tags beginning with "_" are internal, not added to tagname list
    if tag[0] == "_":
      return False;
    # add to list
    self.tagnames.append(tag);
    self.tagnames.sort();
    # add to properties
    self._tagprops[tag] = Source.Property("tag: %s"%tag,
                                                                  lambda src,tag=tag:getattr(src,tag,None) not in [None,False],
                                                                  sources=self.sources,
                                                                  style=self.plotstyles.setdefault('tag:%s'%tag,PlotStyles.DefaultPlotStyle));
    # reform property list
    self._remakePropList();
    return True;
  
  def updateSelection (self):
    """recomputes number of selected sources""";
    self.selprop.computeTotal(self.sources);
    return self.selprop.total;
    
  def setFieldCenter (self,ra0,dec0):
    self.ra0 = ra0;
    self.dec0 = dec0;
    
  def setPrimaryBeam (self,pbexp):
    self.pbexp = pbexp;
    # eval("lambda r,fq:"+self.beam_expr);
  
  def computeSourceLMs (self,proj='sin'):
    """ensures that source lm_sin or lm_ncp coordinates are properly computed and stored in the source items.
    Returns max abs lm offset.""";
    # nothing to do if no source list
    if not self.sources:
      return;
    # first, ensure that field center is set. If not, use model's center of gravity
    update = False;
    if self.ra0 is None:
      self.ra0 = reduce(lambda x,y:x+y,[ src.pos.ra for src in self.sources ])/len(self.sources);
      update = self.updated = True;
    if self.dec0 is None:
      self.dec0 = reduce(lambda x,y:x+y,[ src.pos.dec for src in self.sources ])/len(self.sources);
      update = self.updated = True;
    # now, go through all sources and compute their LMs
    if proj == 'sin':
      coordinate = '_lm_sin';
    elif proj == 'ncp':
      coordinate = '_lm_ncp';
    lm_max = 0.;
    for src in self.sources:
      # for each source, if field center was updated, or if lm_sin or lm_ncp attribute is missing,
      # use the corresponding lm_sin() or lm_ncp() method in its position object to compute lm's.
      lm = getattr(src,coordinate,None);
      if update or not lm:
        lm = getattr(src.pos,coordinate[1:])(self.ra0,self.dec0)
        src.setAttribute(coordinate,lm);
        self.updated = True;
      lm_max = max(lm_max,abs(lm[0]),abs(lm[1]));
    return lm_max;
    
SkyModel.registerClass();
