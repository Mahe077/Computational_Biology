from Bio import Entrez

Entrez.email = "maheshlakshan766@gmail.com"

handle = Entrez.esearch(db="pubmed", term="Alzheimer's")
record = Entrez.read(handle)
pubmed_ids = record["IdList"]

if pubmed_ids:
    handle = Entrez.efetch(db="pubmed", id=pubmed_ids, rettype="xml")
    records = Entrez.read(handle)

    print("Total number of articles:", record['Count'])

    for record in records["PubmedArticle"]:
        authors = ", ".join(author["LastName"] + " " + author["Initials"] for author in record["MedlineCitation"]["Article"]["AuthorList"])
        print(f"Authors: {authors}")
else:
    print("No articles found.")
