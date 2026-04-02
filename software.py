from abc import ABC, abstractmethod
import pandas as pd
import numpy as np
from flask import Flask, render_template, request
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ===== ABSTRACT PARSER =====
class BaseParser(ABC):

    def __init__(self, filepath):
        self.filepath = filepath
        self.data = None

    @abstractmethod
    def parse(self):
        pass

    def get_data(self):
        return self.data


# ===== GO PARSER =====
class GOParser(BaseParser):

    def parse(self):
        terms = []
        current = {}

        with open(self.filepath, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()

                if line == "[Term]":
                    if current:
                        terms.append(current)
                    current = {}

                elif ":" in line:
                    key, value = line.split(":", 1)

                    if key.strip() == "is_a":
                        current.setdefault("is_a", []).append(value.strip())
                    else:
                        current[key.strip()] = value.strip()

            if current:
                terms.append(current)

        self.data = pd.DataFrame(terms)
        return self.data


# ===== GAF PARSER =====
class GAFParser(BaseParser):

    def parse(self):
        rows = []

        with open(self.filepath, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith("!"):
                    continue
                rows.append(line.strip().split("\t"))

        columns = [
            "DB","DB_Object_ID","DB_Object_Symbol","Qualifier",
            "GO_ID","DB_Reference","Evidence_Code","With_From",
            "Aspect","DB_Object_Name","Synonym","DB_Object_Type",
            "Taxon","Date","Assigned_By","Annotation_Extension",
            "Gene_Product_Form_ID"
        ]

        self.data = pd.DataFrame(rows, columns=columns)
        return self.data


# ===== GO TERM =====
class GOTerm:

    def __init__(self, go_id, name, locus):
        self.go_id = go_id
        self.name = name
        self.locus = locus
        self.parents = set()
        self.children = set()

    def add_parent(self, parent):
        self.parents.add(parent)

    def add_child(self, child):
        self.children.add(child)


# ===== ONTOLOGY =====
class Ontology:

    def __init__(self):
        self.terms = {}

    def add_term(self, term):
        self.terms[term.go_id] = term

    def get_term(self, go_id):
        return self.terms.get(go_id)

    def add_relation(self, child_id, parent_id):
        child = self.get_term(child_id)
        parent = self.get_term(parent_id)

        if child and parent:
            child.add_parent(parent)
            parent.add_child(child)

    # ===== NAVIGATION =====

    def get_ancestors(self, go_id):
        term = self.get_term(go_id)
        visited = set()

        def dfs(t):
            for parent in t.parents:
                if parent.go_id not in visited:
                    visited.add(parent.go_id)
                    dfs(parent)

        if term:
            dfs(term)

        return visited

    def get_descendants(self, go_id):
        term = self.get_term(go_id)
        visited = set()

        def dfs(t):
            for child in t.children:
                if child.go_id not in visited:
                    visited.add(child.go_id)
                    dfs(child)

        if term:
            dfs(term)

        return visited

    def get_depth(self, go_id):
        term = self.get_term(go_id)

        if not term or not term.parents:
            return 0

        return 1 + max(self.get_depth(p.go_id) for p in term.parents)


# ===== BUILD ONTOLOGY =====
def build_ontology(go_df):

    ontology = Ontology()

    for _, row in go_df.iterrows():
        term = GOTerm(row['id'], row.get('name'), row.get('locus'))
        ontology.add_term(term)

    for _, row in go_df.iterrows():
        if 'is_a' in row and isinstance(row['is_a'], list):
            for parent in row['is_a']:
                parent_id = parent.split()[0]
                ontology.add_relation(row['id'], parent_id)

    return ontology


# ===== ANNOTATION =====
class Annotation:

    def __init__(self, gene, go_id, evidence):
        self.gene = gene
        self.go_id = go_id
        self.evidence = evidence


# ===== ANNOTATION SET =====
class AnnotationSet:

    def __init__(self):
        self.annotations = []
        self.by_gene = {}

    def add_annotation(self, annotation):

        existing = self.by_gene.get(annotation.gene, [])

        if any(a.go_id == annotation.go_id for a in existing):
            return

        self.annotations.append(annotation)
        self.by_gene.setdefault(annotation.gene, []).append(annotation)

    def get_annotations_for_gene(self, gene):
        return self.by_gene.get(gene, [])


# ===== BUILD ANNOTATIONS =====
def build_annotations(gaf_df):

    aset = AnnotationSet()

    for _, row in gaf_df.iterrows():
        ann = Annotation(
            gene=row['DB_Object_Symbol'],
            go_id=row['GO_ID'],
            evidence=row['Evidence_Code']
        )
        aset.add_annotation(ann)

    return aset


# ===== SIMILARITY =====
class SimilarityMeasure(ABC):

    @abstractmethod
    def compute(self, set1, set2):
        pass


class JaccardSimilarity(SimilarityMeasure):

    def compute(self, set1, set2):
        intersection = len(set1 & set2)
        union = len(set1 | set2)
        return intersection / union if union != 0 else 0


class OverlapSimilarity(SimilarityMeasure):

    def compute(self, set1, set2):
        intersection = len(set1 & set2)
        minimum = min(len(set1), len(set2))
        return intersection / minimum if minimum != 0 else 0


# ===== GENE TERMS =====
def gene_terms(annotation_set, ontology, gene):

    anns = annotation_set.get_annotations_for_gene(gene)
    terms = set()

    for a in anns:
        term = ontology.get_term(a.go_id)
        if term is None:
            continue

        terms.add(a.go_id)
        terms.update(ontology.get_ancestors(a.go_id))

    return terms


def gene_similarity(gene1, gene2, annotation_set, ontology, measure):

    t1 = gene_terms(annotation_set, ontology, gene1)
    t2 = gene_terms(annotation_set, ontology, gene2)

    return measure.compute(t1, t2)


def similarity_matrix(genes, annotation_set, ontology, measure):

    n = len(genes)
    matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            matrix[i, j] = gene_similarity(
                genes[i],
                genes[j],
                annotation_set,
                ontology,
                measure
            )

    return matrix


# ===== GRAPH =====
def plot_go_graph(ontology, go_id):

    G = nx.DiGraph()

    ancestors = ontology.get_ancestors(go_id)
    descendants = ontology.get_descendants(go_id)

    nodes = set(ancestors) | set(descendants)
    nodes.add(go_id)

    for n in nodes:
        term = ontology.get_term(n)
        if term:
            for parent in term.parents:
                if parent.go_id in nodes:
                    G.add_edge(parent.go_id, n)

    colors = ["red" if n == go_id else "lightblue" for n in G.nodes()]

    plt.figure(figsize=(8, 6))
    pos = nx.spring_layout(G)

    nx.draw(G, pos,
            with_labels=True,
            node_color=colors,
            node_size=500,
            font_size=6)

    filepath = "static/graph.png"
    plt.savefig(filepath)
    plt.close()

    return filepath


# ===== FLASK APP =====
app = Flask(__name__)

# LOAD DATA
go_parser = GOParser("go-basic.obo")
go_df = go_parser.parse()

gaf_parser = GAFParser("goa_human.gaf")
gaf_df = gaf_parser.parse()

ontology = build_ontology(go_df)
annotations = build_annotations(gaf_df)


@app.route("/")
def home():
    return render_template(
        "home.html",
        genes=list(annotations.by_gene.keys()),
        gos=list(ontology.terms.keys())
    )


@app.route("/gene")
def gene_page():

    gene = request.args.get("gene")
    anns = annotations.get_annotations_for_gene(gene)

    results = []

    for ann in anns[:20]:
        term = ontology.get_term(ann.go_id)
        if term is None:
            continue

        ancestors = list(ontology.get_ancestors(ann.go_id))[:5]

        results.append({
            "go_id": ann.go_id,
            "name": term.name,
            "ancestors": ancestors
        })

    return render_template("gene.html", gene=gene, results=results)


@app.route("/go")
def go_page():

    go_id = request.args.get("go")
    term = ontology.get_term(go_id)

    if term is None:
        return render_template("go.html", go_id=go_id, error=True)

    graph_path = plot_go_graph(ontology, go_id)

    return render_template(
        "go.html",
        go_id=go_id,
        name=term.name,
        ancestors=list(ontology.get_ancestors(go_id))[:10],
        descendants=list(ontology.get_descendants(go_id))[:10],
        depth=ontology.get_depth(go_id),
        graph=graph_path,
        error=False
    )


@app.route("/similarity")
def similarity_page():

    gene1 = request.args.get("gene1")
    gene2 = request.args.get("gene2")
    measure_name = request.args.get("measure")

    measure = OverlapSimilarity() if measure_name == "overlap" else JaccardSimilarity()

    sim = gene_similarity(gene1, gene2, annotations, ontology, measure)
    matrix = similarity_matrix([gene1, gene2], annotations, ontology, measure)

    return render_template(
        "similarity.html",
        gene1=gene1,
        gene2=gene2,
        similarity=sim,
        matrix=matrix,
        measure=measure_name,
        genes=list(annotations.by_gene.keys())
    )


if __name__ == "__main__":
    app.run(debug=True)
