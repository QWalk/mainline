from __future__ import print_function
import yaml




def print_entry_md(f,en):
    print("*******************************",file=f)
    print("*******************************",file=f)
    
    print("##"+en['title'],file=f)
    print("Keyword:",en['keyword'],file=f)
    print("\n\n",file=f)
    
    print("*Description:*",en['description'],'\n',file=f)
    
    print("*Required keywords*\n\n",file=f)
    if len(en['required']) > 0:
      print("Keyword|Type|Description",file=f)
      print("------------ | ------------- | ------------",file=f)
      for i in en['required']:
         kw=i['keyword']
         print(kw,'|',i['type'],'|',i['description'].replace('\n',' '),file=f)
    else: 
      print("None",file=f)
    print("\n",file=f)

    print("*Optional keywords*\n\n",file=f)
    if len(en['optional']) > 0:

      print("Keyword|Type|Default|Description",file=f)
      print("------------ | ------------- | ------------- | ------------",file=f)
      for i in en['optional']:
         kw=i['keyword']
         print(kw,'|',i['type'],'|',i['default'],'|',i['description'].replace('\n',' '),file=f)
    else:
      print("None",file=f)
    print("\n\n",file=f)
    if 'advanced' in en.keys():
      print("*Advanced keywords*\n\n",file=f)
      print("Keyword|Type|Default|Description",file=f)
      print("------------ | ------------- | ------------- | ------------",file=f)
      for i in en['advanced']:
         kw=i['keyword']
         print(kw,'|',i['type'],'|',i['default'],'|',i['description'].replace('\n',' '),file=f)

    print("\n",file=f)


def print_index(f):
  print("""# QWalk documentation
The main components of a QMC (or just quantum!) calculation are:

1. The [Hamiltonian](Hamiltonian). This is the contents of the SYSTEM section.
2. The [Trial Wavefunction](Wavefunction). This is the contents of the TRIALFUNC section.
3. The [Method](Method) used to evaluate properties. This is the contents of the METHOD section
4. The [Observables](Average generator) calculated. These are generally inserted in the METHOD section with the keyword AVERAGE.

The input of QWalk is set up to provide those four components of a calculation in a modular way. To assist with this, there are several important helper components:

1. Simple [Basis function](Basis function)s to represent 3-dimensional functions in real space.
2. [One-particle orbitals](Orbital) made up of linear combinations of basis functions.
3. [Converter](Converter) programs to translate wave functions from DFT and quantum chemistry codes to QWalk format.

## Input rules

* Whitespace is not important
* You can comment a line with '#'
* Keywords are case-insensitive, but string arguments are case-sensitive
* There are three types of keywords
     * flags: no argument. For example, 'TMOVES'.
     * floats, integers and strings: a single argument. For example, 'VMC_NCONFIG 20000'
     * sections: a collection of keywords denoted by curly braces { }. For example, 'METHOD { VMC NBLOCK 50 NDECORR 3 AVERAGE { SK } }'
* The order of keywords is not important unless a keyword is repeated. For example, multiple METHOD sections in an input file will execute them in order.

""", file=f)


if __name__=="__main__":
  import glob
  files=glob.glob("*.yml")
  entries=[]
  for fi in files:
    entries.append(yaml.load(open(fi)))
  
  #TeX-AMS-MML_HTMLorMML
  print_index(open("../mkdocs/docs/index.md",'w'))

  f=open("../mkdocs/mkdocs.yml",'w')
  f.write("site_name: QWalk documentation\n")
  f.write("theme: readthedocs\n")
  f.write("""extra_javascript: 
    - http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML

markdown_extensions:
    - mdx_math
\n""")
  
  f.write("pages:\n")
  f.write(" - About: index.md\n")
  categories=[]
  for en in entries:
    if en['type']=="Category":
      categories.append(en['name'])
      #f.write(" - "+en['name']+": "+en['name']+".md")
  #f.close()

  for c in categories:
    f1=open("../mkdocs/docs/"+c+".md",'w')
    f.write(" - "+c+": "+c+".md\n")
    #f.write(" - "+c+":\n")
    
    for en in entries:
      if en['type']=='Entry' and en['is_a']==c:
        #f.write("    - "+en['name']+": "+en['name']+".md\n")
        #f1=open("../mkdocs/docs/"+en['name']+".md",'w')
        print_entry_md(f1,en)
      

  

