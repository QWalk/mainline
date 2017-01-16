
all: 
	cd database; python gen_md.py
	cd mkdocs; mkdocs build --clean

tgz: all
	cd mkdocs/site; tar -czf ../../docs.tgz * 

serve: all
	cd mkdocs; mkdocs serve
