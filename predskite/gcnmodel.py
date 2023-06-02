import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv


class GCN(torch.nn.Module):
    """Graph Convolutional Network (GCN) for node classification."""
    
    def __init__(self, input_dim=1, hidden_dim=16, output_dim=1):
        """Initialize GCN with two GCNConv layers."""
        super(GCN, self).__init__()
        self.conv1 = GCNConv(input_dim, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, output_dim)

    def forward(self, data):
        """Define forward pass for GCN."""
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, training=self.training)
        x = self.conv2(x, edge_index)
        return torch.mean(x).view(1)


def train_model(model, train_data, optimizer):
    """Train the model."""
    model.train()
    for epoch in range(200):
        for data in train_data:
            data = data.to(device)
            optimizer.zero_grad()
            output = model(data)
            loss = F.mse_loss(output, data.y)
            loss.backward()
            optimizer.step()


def evaluate_model(model, test_data):
    """Evaluate the model and return the test loss."""
    model.eval()
    test_loss = 0
    for data in test_data:
        data = data.to(device)
        with torch.no_grad():
            output = model(data)
        test_loss += F.mse_loss(output, data.y).item()
    return test_loss
