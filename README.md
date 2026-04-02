# Gene Ontology Web Application

## Description

This project is a web application built with Flask to explore the Gene Ontology and human gene annotations.

It allows users to:

* view gene annotations
* explore GO terms (ancestors, descendants, depth)
* compute similarity between genes
* visualize a small GO graph

---

## Project Structure

```
project/
│
├── software.py
├── templates/
│   ├── home.html
│   ├── gene.html
│   ├── go.html
│   └── similarity.html
│
├── static/
│   └── graph.png
│
└── .gitignore
```

---

## Setup

Download the following files:

* [GO ontology (OBO)](https://current.geneontology.org/ontology/go-basic.obo)
* [Human gene annotations (GAF)](https://current.geneontology.org/annotations/goa_human.gaf.gz)

Place them in the same folder as `software.py`.
___

## File Placement

After downloading the `.obo` and `.gaf` files, place them in the same directory as `software.py`.

Do not move or rename the `templates/` and `static/` folders, as they are required by Flask.

---

## Requirements

Install dependencies:

```
pip install flask pandas numpy networkx matplotlib
```

---

## Run the Application

```
python software.py
```

Open the browser at:

```
http://127.0.0.1:5000
```

---

## Pages

* **Home (`home.html`)**
  Main page to select genes and GO terms

* **Gene (`gene.html`)**
  Shows annotations for a selected gene

* **GO (`go.html`)**
  Shows details of a GO term and its graph

* **Similarity (`similarity.html`)**
  Computes similarity between two genes

---

## Notes

* The `.obo` and `.gaf` files are not included in the repository
* They must be downloaded manually before running the application
* The graph image is generated automatically in the `static/` folder

---

## Author

Gianmarco Macrì (https://github.com/GianmarcoMacri)
Ahmad Shakaroun (https://github.com/Ahmadshakaroun)
