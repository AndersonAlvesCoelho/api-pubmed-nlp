from pipeline.data_fetch import fetch_pubmed_articles
from pipeline.preprocessing import preprocess_abstracts
from pipeline.summarizer import summarize_abstracts, summarize_general
from pipeline.summarizer_bert import summarize_abstracts_bert, summarize_general_bert

def analyze_articles(pubmed_ids: list[str], email: str, use_bert: bool = True) -> dict:
    """
    Pipeline principal — busca artigos, processa textos e retorna resumos.
    """
    # 1. Buscar artigos
    articles = fetch_pubmed_articles(pubmed_ids, email)
    if not articles:
        return {"error": "Nenhum artigo encontrado."}

    # 2. Pré-processamento
    processed = preprocess_abstracts(articles)

    # 3. Resumo individual e geral
    if use_bert:
        summarized = summarize_abstracts_bert(processed)
        general_summary = summarize_general_bert([a["processed"] for a in summarized])
    else:
        summarized = summarize_abstracts(processed)
        general_summary = summarize_general([a["abstract"] for a in summarized])

    return {
        "articles": summarized,
        "general_summary": general_summary,
        "model": "BERT" if use_bert else "LSA"
    }
