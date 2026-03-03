"""
Tests for network analysis functions (requires networkx).
"""

import pytest

from llps.network import _classify_network_nodes, _compute_network_metrics


def test_classify_network_nodes() -> None:
    pytest.importorskip("networkx")
    import networkx as nx
    G = nx.Graph()
    G.add_node('A', pLLPS=0.9)
    G.add_node('B', pLLPS=0.8)
    G.add_node('C', pLLPS=0.3)
    G.add_node('D', pLLPS=None)
    high, low, unknown = _classify_network_nodes(G, high_threshold=0.7)
    assert 'A' in high
    assert 'B' in high
    assert 'C' in low
    assert 'D' in unknown


def test_compute_network_metrics_empty_graph() -> None:
    pytest.importorskip("networkx")
    import networkx as nx
    G = nx.Graph()
    results = _compute_network_metrics(G, [], [], [])
    assert results['total_nodes'] == 0
    assert results['total_edges'] == 0
    assert results['enrichment_ratio'] == 0.0


def test_compute_network_metrics_simple() -> None:
    pytest.importorskip("networkx")
    import networkx as nx
    G = nx.Graph()
    G.add_node('A', pLLPS=0.9)
    G.add_node('B', pLLPS=0.8)
    G.add_node('C', pLLPS=0.3)
    G.add_edge('A', 'B')
    G.add_edge('A', 'C')
    results = _compute_network_metrics(G, ['A', 'B'], ['C'], [])
    assert results['total_nodes'] == 3
    assert results['total_edges'] == 2
    assert results['high_high_interactions'] == 1
    assert results['high_low_interactions'] == 1
    assert results['low_low_interactions'] == 0
