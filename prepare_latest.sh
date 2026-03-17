#!/bin/bash

# 1. Encontrar a pasta mais recente
LATEST_DIR=$(find biocypher-out -mindepth 1 -maxdepth 1 -type d ! -name "latest" | sort -r | head -n 1)

if [ -z "$LATEST_DIR" ]; then
    echo "❌ Nenhuma pasta encontrada em biocypher-out!"
    exit 1
fi

echo "📂 Pasta detetada: $LATEST_DIR"

# 2. Atualizar o .env para o Docker
echo "IMPORT_DIR=$LATEST_DIR" > .env
echo "⚙️ .env atualizado."

# 3. Criar o ficheiro .sh DO ZERO
SCRIPT_PATH="$LATEST_DIR/neo4j-admin-import-call.sh"

# --- A MUDANÇA ESTÁ AQUI: Forçar a remoção se o ficheiro já existir ---
rm -f "$SCRIPT_PATH" 

# Criar o ficheiro de novo (usando > para criar/limpar)
echo "#!/bin/bash" > "$SCRIPT_PATH"
# ---------------------------------------------------------------------

{
    echo "neo4j-admin import \\"
    echo "  --database=neo4j \\"
    echo "  --delimiter=\"\t\" \\"
    echo "  --array-delimiter=\"|\" \\"
    echo "  --quote=\"'\" \\"
    echo "  --force=true \\"
    echo "  --skip-bad-relationships=true \\"
    echo "  --skip-duplicate-nodes=true \\"

    # Adicionar todos os NÓS
    for header in "$LATEST_DIR"/*-header.csv; do
        filename=$(basename "$header")
        base=${filename%-header.csv}
        if [[ ! $base =~ (Has|In|ReactionMetabolite|GeneReaction) ]]; then
            echo "  --nodes=\"/data/build2neo/$base-header.csv,/data/build2neo/$base-part.*\" \\"
        fi
    done

    # Adicionar todas as ARESTAS
    for header in "$LATEST_DIR"/*-header.csv; do
        filename=$(basename "$header")
        base=${filename%-header.csv}
        if [[ $base =~ (Has|In|ReactionMetabolite|GeneReaction) ]]; then
            echo "  --relationships=\"/data/build2neo/$base-header.csv,/data/build2neo/$base-part.*\" \\"
        fi
    done
} >> "$SCRIPT_PATH" # Escreve tudo de uma vez para evitar múltiplos acessos

# Remover a última barra invertida
sed -i '$ s/ \\//' "$SCRIPT_PATH"
# Converter CRLF para LF (importante no Windows!)
sed -i 's/\r$//' "$SCRIPT_PATH"

chmod +x "$SCRIPT_PATH"
echo "✨ Ficheiro $SCRIPT_PATH gerado com sucesso!"