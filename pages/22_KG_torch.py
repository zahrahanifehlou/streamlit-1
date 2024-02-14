from torch.optim import Adam
from torch import cuda
from torchkge.models import TransEModel
from torchkge.sampling import BernoulliNegativeSampler
from torchkge.utils import MarginLoss, DataLoader
from torchkge.utils.datasets import load_fb15k
from torchkge.evaluation import LinkPredictionEvaluator
from torchkge.utils.pretrained_models import load_pretrained_transe

from tqdm.autonotebook import tqdm

import streamlit as st

emb_dim = 100
lr = 0.0004
margin = 0.5
n_epochs = 1000
b_size = 32768

# Load dataset
kg_train, kg_val, kg_test = load_fb15k()
st.write(kg_train.get_df())
st.write(kg_test.get_df())
st.write(kg_val.get_df())
# Define the model and criterion
# Define the model and criterion
# model = TransEModel(emb_dim, kg_train.n_ent, kg_train.n_rel, dissimilarity_type="L2")
# criterion = MarginLoss(margin)
model = load_pretrained_transe("fb15k", 100)
# Move everything to CUDA if available
if cuda.is_available():
    cuda.empty_cache()
    model.cuda()
#     criterion.cuda()

# # Define the torch optimizer to be used
# optimizer = Adam(model.parameters(), lr=lr, weight_decay=1e-5)

# sampler = BernoulliNegativeSampler(kg_train)
# dataloader = DataLoader(kg_train, batch_size=b_size, use_cuda="all")

# iterator = tqdm(range(n_epochs), unit="epoch")
# for epoch in iterator:
#     running_loss = 0.0
#     for i, batch in enumerate(dataloader):
#         h, t, r = batch[0], batch[1], batch[2]
#         n_h, n_t = sampler.corrupt_batch(h, t, r)

#         optimizer.zero_grad()

#         # forward + backward + optimize
#         pos, neg = model(h, t, r, n_h, n_t)
#         loss = criterion(pos, neg)
#         loss.backward()
#         optimizer.step()

#         running_loss += loss.item()
#     iterator.set_description(
#         "Epoch {} | mean loss: {:.5f}".format(epoch + 1, running_loss / len(dataloader))
#     )
#     st.write("epoch", epoch + 1)
#     st.write("running loss", running_loss / len(dataloader))

# model.normalize_parameters()


evaluator = LinkPredictionEvaluator(model, kg_test)
evaluator.evaluate(b_size=32)
st.write(evaluator.print_results())
