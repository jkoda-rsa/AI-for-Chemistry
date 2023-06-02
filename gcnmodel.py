import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv


class GCN(torch.nn.Module):    
    def __init__(self, input_dim=1, hidden_dim=16, output_dim=1):
        super(GCN, self).__init__()
        self.conv1 = GCNConv(input_dim, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, output_dim)

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, training=self.training)
        x = self.conv2(x, edge_index)
        return torch.mean(x).view(1)


def train_model(model, train_data, optimizer):
    model.train()
    for epoch in range(20):
        for data in train_data:
            data = data.to(device)
            optimizer.zero_grad()
            output = model(data)
            loss = F.mse_loss(output, data.y)
            train_loss = loss
            loss.backward()
            optimizer.step()
    print('Training Loss:', loss)
            

def evaluate_model(model, test_data):
    model.eval()
    test_loss = 0
    for data in test_data:
        data = data
        with torch.no_grad():
            output = model(data)
        test_loss += F.mse_loss(output, data.y).item()
    return test_loss

if __name__ == "__main__":
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    model = GCN().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)
    
    train_model(model, train_data, optimizer)
    
    test_loss = evaluate_model(model, test_data)
    print(f'Test Loss: {test_loss}')
