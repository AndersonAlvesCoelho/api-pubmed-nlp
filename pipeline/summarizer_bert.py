from transformers import pipeline

# Carregamos um modelo de sumarização baseado em BERT ou T5
# (T5 e BART são arquiteturas encoder-decoder baseadas em BERT)
summarizer = pipeline("summarization", model="facebook/bart-large-cnn")

def summarize_with_bert(text: str, max_length: int = 200, min_length: int = 60) -> str:
    """
    Gera um resumo semântico usando modelo BERT-like.
    """
    if not text or len(text.split()) < 30:
        return text  # evita erro com textos muito curtos

    summary = summarizer(
        text,
        max_length=max_length,
        min_length=min_length,
        do_sample=False
    )
    return summary[0]['summary_text']

def summarize_abstracts_bert(articles: list[dict]) -> list[dict]:
    """
    Gera resumo para cada artigo individualmente com BERT.
    """
    for a in articles:
        a["summary"] = summarize_with_bert(a["processed"])
    return articles

def summarize_general_bert(abstracts: list[str]) -> str:
    """
    Cria um resumo geral combinando todos os abstracts.
    """
    combined = " ".join(abstracts)
    return summarize_with_bert(combined, max_length=400, min_length=120)
