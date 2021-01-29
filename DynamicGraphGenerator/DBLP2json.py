#!/usr/bin/python
## wget -N http://dblp.uni-trier.de/xml/dblp.xml.gz
## then run this script
import codecs, collections, json, gzip, os, sys, xml.sax

xml_filename = 'data/dblp.xml.gz'
json_gz_filename = 'data/dblp.conf_journal.json.gz'
tmp_filename = 'data/tmp_dblp.conf_journal.json.gz'
report_frequency = 10000
# areas = ['theory'] # Choice for the area:


class DBLPHandler (xml.sax.ContentHandler):
    # removing "www" from the list below as www corresponds to the homepage of an author
    # removing "phdthesis"  and "mastersthesis"
    # Removing proceedings, which collects the proceeding of a conference
    # Removing 'incollection' as we only want to focus on articles and conference versions
    # papertypes = set (['article', 'book', 'inproceedings', 'incollection', 'www', 'proceedings', 'phdthesis', 'mastersthesis'])
    papertypes = {'article', 'inproceedings'}
    # conferences = {}
    # conferences['theory'] = ['/soda/','/stoc/','/focs/']
    def __init__ (self, out):
        self.out = out
        self.paper = None
        self.authors = []
        self.year = None
        self.text = ''
        # self.booktitle = ''
        # self.journal = ''
        self.papercount = 0
        self.edgecount = 0
    def startElement (self, name, attrs):
        if name in self.papertypes:
            self.paper = str (attrs['key'])
            self.authors = []
            self.year = None
            self.journal = ''
            self.booktitle = ''
        elif name in ['author', 'year','booktitle','journal']:
            self.text = ''
    def endElement (self, name):
        if name == 'author':
            self.authors.append (self.text)
        if name == 'year':
            self.year = int (self.text.strip ())
        # if name == 'booktitle':
        #     self.booktitle = self.text
        # if name == 'journal':
        #     self.journal = self.text
        elif name in self.papertypes:
            # if name == 'inproceedings' and self.match_conference():
            #     self.write_paper ()
            # if name == 'article' and self.match_journal():
            self.write_paper ()
            self.paper = None

    def write_paper (self):
        if self.papercount:
            self.out.write (',\n')
        self.papercount += 1
        self.edgecount += len (self.authors)
        json.dump ([self.paper, self.authors, self.year], self.out)
        if self.papercount % report_frequency == 0:
            print ('... processed %d papers, %d edges so far ...' % (self.papercount, self.edgecount))
            sys.stdout.flush ()
    def characters (self, chars):
        self.text += chars

    # def match_conference(self):
    #     for area in areas:
    #         for conf in  self.conferences[area]:
    #             if conf in self.paper:
    #                 print(self.paper,self.year)
    #                 return True
    #     return False





def force ():
    print ('** Parsing XML...')

    xmlfile = gzip.GzipFile (xml_filename, 'r')
    out = gzip.open (tmp_filename, 'wt')
    out.write ('[\n')
    dblp = DBLPHandler (out)
    parser = xml.sax.parse (xmlfile, dblp)
    out.write ('\n]\n')
    out.close ()
    os.rename (tmp_filename, json_gz_filename)

    print ('-- %d papers, %d edges' % (dblp.papercount, dblp.edgecount))

def main (parse_args = False):
    try:
        need = (os.stat (xml_filename).st_mtime >= os.stat (json_gz_filename).st_mtime)
    except OSError:
        need = True
    if parse_args and len (sys.argv) > 1:
        need = True
    if need:
        force ()

def open ():
    main ()
    return gzip.GzipFile (json_gz_filename, 'r')

def papers ():
    for line in open ():
        if line.strip () in '[]': continue
        line = line.rstrip ().rstrip (',')
        yield json.loads (line)

if __name__ == '__main__': main (True)

# Security: IEEE Security and Privacy symposium; USENIX Security; ACM Computer and Communication
# Security; Network and Distributed systems security symposium (NDSS); Annual Computer Security
# Application conference

# Systems: Operating System Design and Implementation; Symposium on OS Principles;
# NSDI (Networked Systems Design and Implementations); FAST File and Storage Technologies

#