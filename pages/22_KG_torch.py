from torch.optim import Adam
from torch import cuda
from torchkge.models import TransEModel
from torchkge.sampling import BernoulliNegativeSampler
from torchkge.utils import MarginLoss, DataLoader
from torchkge.utils.datasets import load_fb15k
from torchkge.evaluation import LinkPredictionEvaluator
from torchkge.utils.pretrained_models import load_pretrained_transe
import pandas as pd
from torchkge.data_structures import KnowledgeGraph
from torchkge.inference import RelationInference
from tqdm.autonotebook import tqdm
import plotly.express as px
import streamlit as st
from sklearn.model_selection import train_test_split

lr = 0.0004
margin = 0.5
n_epochs = 1500
b_size = 1000


@st.cache_data
def load_data():
    primekg = pd.read_csv("../kg.csv", low_memory=False)
    return primekg


# Load dataset
kg = load_data()
kg_train = kg[["relation", "x_name", "y_name"]]
# kg_train2, kg_val, kg_test = load_fb15k()


# st.write(kg_train2.get_df())
kg_train = kg_train.sample(100000)
# df.rename(columns={"A": "a", "B": "c"})
kg_train = kg_train.rename(
    columns={"x_name": "from", "y_name": "to", "relation": "rel"}
)
# kg_train, kg_test = train_test_split(kg_train, test_size=0.3, stratify=kg_train["rel"])
st.write("size kg", (len(kg_train["from"].unique()), len(kg_train.rel.unique())))
st.write(kg_train.sample(20))
kg_t = KnowledgeGraph(kg_train)
st.warning("DONE CONVERSION!")
emb_dim = len(kg_train.rel.unique())
emb_dim = 5
# exit(0)
# st.write(kg_val.get_df())
# Define the model and criterion
# Define the model and criterion
model = TransEModel(emb_dim, kg_t.n_ent, kg_t.n_rel, dissimilarity_type="L2")
criterion = MarginLoss(margin)
# model = load_pretrained_transe("fb15k", 100)
# Move everything to CUDA if available
if cuda.is_available():
    cuda.empty_cache()
    model.cuda()
    criterion.cuda()

# Define the torch optimizer to be used
optimizer = Adam(model.parameters(), lr=lr, weight_decay=1e-5)

sampler = BernoulliNegativeSampler(kg_t)
dataloader = DataLoader(kg_t, batch_size=b_size, use_cuda="all")

losses = []
iterator = tqdm(range(n_epochs), unit="epoch")
for epoch in iterator:
    running_loss = 0.0
    for i, batch in enumerate(dataloader):
        h, t, r = batch[0], batch[1], batch[2]
        n_h, n_t = sampler.corrupt_batch(h, t, r)

        optimizer.zero_grad()

        # forward + backward + optimize
        pos, neg = model(h, t, r, n_h, n_t)
        loss = criterion(pos, neg)
        loss.backward()
        optimizer.step()

        running_loss += loss.item()
    iterator.set_description(
        "Epoch {} | mean loss: {:.5f}".format(epoch + 1, running_loss / len(dataloader))
    )
    losses.append(running_loss / len(dataloader))
    model.normalize_parameters()
    # st.write("epoch", epoch + 1)
    # st.write("running loss", running_loss / len(dataloader))
model.normalize_parameters()
# ent_emb, rel_emb = model.get_embeddings()

# for a in emb:
#     st.write(a.shape)
df_loss = pd.DataFrame()
df_loss["loss"] = losses
# st.write(ent_emb)
st.plotly_chart(px.line(df_loss, x=df_loss.index, y="loss"))
# model.inference_prepare_candidates()
# a = RelationInference(model, ent_emb, ent_emb.indices())
# model.RelationInference.evaluate(b_size=b_size)
evaluator = LinkPredictionEvaluator(model, kg_t)
evaluator.evaluate(b_size=32)
st.write(evaluator.print_results())
