import spextra as sp
from spextra import Database

db = Database()
print(db.liblist)

print(sp.database.SpecLibrary(db.liblist[0]).template_names)