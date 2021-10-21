make ssa_exp
category="$1"
echo "Dataset"
echo "$category"
echo "TopK from config.txt"
topK="$2"
echo "$topK"

configK="$topK \/\/ C"
sed -i'' -e "1s/.*/$configK/" config.txt

./bin/ssa -l -c "$category"

echo "[HNSW] Search Start"
cd ip-nsw
bash run.sh "$category" "$topK"
echo "[HNSW] Search End"

cd ..
mkdir -p figure
mv "$category"_*"$topK".txt figure/
