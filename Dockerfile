# Imagem base
FROM python:3.12-slim

# Define o diretório de trabalho dentro do container
WORKDIR /app

# Evita prompts interativos e instala dependências do sistema necessárias
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    libssl-dev \
    libffi-dev \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

# Copia o arquivo de dependências e instala pacotes Python
COPY requirements.txt . 
RUN pip install --no-cache-dir --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

# Instala spaCy e baixa o modelo en_core_web_sm
RUN pip install --no-cache-dir spacy==3.8.7
RUN python -m spacy download en_core_web_sm

# Copia todo o código da API para o container
COPY . .

# Expõe a porta que a API vai rodar
EXPOSE 8000

# Comando para rodar o uvicorn
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000", "--reload"]
