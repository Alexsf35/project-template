#!/bin/bash
IMPORT_FILE="/data/build2neo/neo4j-admin-import-call.sh"

if [ -f "$IMPORT_FILE" ]; then
    echo "🧹 A traduzir o script para o formato Linux..."
    
    # 1. Converte CRLF (Windows) para LF (Linux)
    sed -i 's/\r$//' "$IMPORT_FILE"
    
    # 2. A MAGIA: Remove caminhos absolutos do Windows (C:\... ou c:/...)
    # Esta linha corta tudo desde o "C:" até à última barra antes do nome do ficheiro CSV
    sed -i 's|[a-zA-Z]:[\\/][^",]*[\\/]||g' "$IMPORT_FILE"
    
    # 3. Corrige o executável: usa o comando global
    sed -i 's|bin/neo4j-admin|neo4j-admin|g' "$IMPORT_FILE"
    
    # 4. Corrige o delimitador para TAB
    sed -i 's/--delimiter="\\t"/--delimiter="TAB"/g' "$IMPORT_FILE"
    sed -i 's/--delimiter="\/t"/--delimiter="TAB"/g' "$IMPORT_FILE"

    chmod +x "$IMPORT_FILE"
    
    # 5. Entra na pasta dos dados e executa
    cd /data/build2neo
    echo "🚀 A iniciar a importação em $(pwd)..."
    bash "$IMPORT_FILE"
else
    echo "⚠️ Script de importação não encontrado!"
fi

# Sinaliza ao Docker que terminámos
neo4j start
sleep 10
neo4j stop
