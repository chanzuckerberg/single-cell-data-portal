const DI_AND_MESENCEPHALON = "Di- and mesencephalon";
const CHOLINERGIC_MONOAMINERGIC = "Cholinergic, monoaminergic...";
const IMMATURE_NEURAL = "Immature neural";
const NEURAL_CREST_LIKE_GLIA = "Neural crest-like glia";
const PERIPHERAL_SENSORY = "Peripheral sensory";
const SPINAL_CORD = "Spinal cord";
const TELENCEPHALON_INTERNEURONS = "Telencephalon interneurons";
const TELENCEPHALON_PROJECTING = "Telencephalon projecting";

export const TISSUE = [
  { id: "brain", name: "Brain" },
  {
    id: "Astroependymal",
    name: "Astroependymal",
    parentId: "brain",
  },
  {
    id: "Cerebellum",
    name: "Cerebellum",
    parentId: "brain",
  },
  {
    id: CHOLINERGIC_MONOAMINERGIC,
    name: CHOLINERGIC_MONOAMINERGIC,
    parentId: "brain",
  },
  {
    id: DI_AND_MESENCEPHALON,
    name: DI_AND_MESENCEPHALON,
    parentId: "brain",
  },
  {
    id: "Enteric",
    name: "Enteric",
    parentId: "brain",
  },
  {
    id: "Hindbrain",
    name: "Hindbrain",
    parentId: "brain",
  },
  {
    id: IMMATURE_NEURAL,
    name: IMMATURE_NEURAL,
    parentId: "brain",
  },
  {
    id: "Immune",
    name: "Immune",
    parentId: "brain",
  },
  {
    id: NEURAL_CREST_LIKE_GLIA,
    name: NEURAL_CREST_LIKE_GLIA,
    parentId: "brain",
  },
  {
    id: "Oligodendrocytes",
    name: "Oligodendrocytes",
    parentId: "brain",
  },
  {
    id: PERIPHERAL_SENSORY,
    name: PERIPHERAL_SENSORY,
    parentId: "brain",
  },
  {
    id: SPINAL_CORD,
    name: SPINAL_CORD,
    parentId: "brain",
  },
  {
    id: "Sympathetic",
    name: "Sympathetic",
    parentId: "brain",
  },
  {
    id: TELENCEPHALON_INTERNEURONS,
    name: TELENCEPHALON_INTERNEURONS,
    parentId: "brain",
  },
  {
    id: TELENCEPHALON_PROJECTING,
    name: TELENCEPHALON_PROJECTING,
    parentId: "brain",
  },
  {
    id: "Vascular",
    name: "Vascular",
    parentId: "brain",
  },
  // Di- and mesencephalon children
  {
    id: "Excitatory",
    name: "Excitatory",
    parentId: DI_AND_MESENCEPHALON,
  },
  {
    id: "Inhibitory",
    name: "Inhibitory",
    parentId: DI_AND_MESENCEPHALON,
  },
  // Excitatory cells
  {
    id: "CR",
    name: "CR",
    parentId: "Excitatory",
  },
  {
    id: "DECHO2",
    name: "DECHO2",
    parentId: "Excitatory",
  },
  {
    id: "DEGLU1",
    name: "DEGLU1",
    parentId: "Excitatory",
  },
  {
    id: "DEGLU2",
    name: "DEGLU2",
    parentId: "Excitatory",
  },
  {
    id: "DEGLU3",
    name: "DEGLU3",
    parentId: "Excitatory",
  },
  {
    id: "DEGLU4",
    name: "DEGLU4",
    parentId: "Excitatory",
  },
  {
    id: "DEGLU5",
    name: "DEGLU5",
    parentId: "Excitatory",
  },
  {
    id: "HBGLU1",
    name: "HBGLU1",
    parentId: "Excitatory",
  },
  {
    id: "HBGLU2",
    name: "HBGLU2",
    parentId: "Excitatory",
  },
  {
    id: "HBGLU3",
    name: "HBGLU3",
    parentId: "Excitatory",
  },
  {
    id: "MEGLU1",
    name: "MEGLU1",
    parentId: "Excitatory",
  },
  {
    id: "MEGLU2",
    name: "MEGLU2",
    parentId: "Excitatory",
  },
  {
    id: "MEGLU3",
    name: "MEGLU3",
    parentId: "Excitatory",
  },
  {
    id: "MEGLU4",
    name: "MEGLU4",
    parentId: "Excitatory",
  },
  {
    id: "MEGLU5",
    name: "MEGLU5",
    parentId: "Excitatory",
  },
  {
    id: "MEGLU6",
    name: "MEGLU6",
    parentId: "Excitatory",
  },
  {
    id: "MEGLU7",
    name: "MEGLU7",
    parentId: "Excitatory",
  },
  {
    id: "MEGLU8",
    name: "MEGLU8",
    parentId: "Excitatory",
  },
  {
    id: "MEGLU9",
    name: "MEGLU9",
    parentId: "Excitatory",
  },
  {
    id: "MEGLU10",
    name: "MEGLU10",
    parentId: "Excitatory",
  },
  {
    id: "MEGLU11",
    name: "MEGLU11",
    parentId: "Excitatory",
  },
];

// GENES
export const NDNF = {
  data: [
    {
      id: "Astroependymal",
      negEntropy: 0,
    },
    {
      id: "Cerebellum",
      negEntropy: 0,
    },
    {
      id: CHOLINERGIC_MONOAMINERGIC,
      negEntropy: 0,
    },
    {
      id: DI_AND_MESENCEPHALON,
      negEntropy: 0.7,
    },
    {
      id: "Enteric",
      negEntropy: 0,
    },
    {
      id: "Hindbrain",
      negEntropy: 0,
    },
    {
      id: IMMATURE_NEURAL,
      negEntropy: 0,
    },
    {
      id: "Immune",
      negEntropy: 0,
    },
    {
      id: NEURAL_CREST_LIKE_GLIA,
      negEntropy: 0,
    },
    {
      id: "Oligodendrocytes",
      negEntropy: 0,
    },
    {
      id: PERIPHERAL_SENSORY,
      negEntropy: 0,
    },
    {
      id: SPINAL_CORD,
      negEntropy: 0,
    },
    {
      id: "Sympathetic",
      negEntropy: 0,
    },
    {
      id: TELENCEPHALON_INTERNEURONS,
      negEntropy: 0,
    },
    {
      id: TELENCEPHALON_PROJECTING,
      negEntropy: 0,
    },
    {
      id: "Vascular",
      negEntropy: 0,
    },
    // Di- and mesencephalon children
    {
      id: "Excitatory",
      negEntropy: 0.7,
      parentId: DI_AND_MESENCEPHALON,
    },
    {
      id: "Inhibitory",
      negEntropy: 0.7,
      parentId: DI_AND_MESENCEPHALON,
    },
    // Excitatory cells
    {
      id: "CR",
      parentId: "Excitatory",
      proportionalExpression: 0.9,
      relativeExpression: 0.85,
    },
    {
      id: "DECHO2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU4",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU5",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU4",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU5",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU6",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU7",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU8",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU9",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU10",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU11",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
  ],
  name: "Ndnf",
};

export const NHLH2 = {
  data: [
    {
      id: "Astroependymal",
      negEntropy: 0,
    },
    {
      id: "Cerebellum",
      negEntropy: 0,
    },
    {
      id: CHOLINERGIC_MONOAMINERGIC,
      negEntropy: 0,
    },
    {
      id: DI_AND_MESENCEPHALON,
      negEntropy: 0.7,
    },
    {
      id: "Enteric",
      negEntropy: 0,
    },
    {
      id: "Hindbrain",
      negEntropy: 0,
    },
    {
      id: IMMATURE_NEURAL,
      negEntropy: 0,
    },
    {
      id: "Immune",
      negEntropy: 0,
    },
    {
      id: NEURAL_CREST_LIKE_GLIA,
      negEntropy: 0,
    },
    {
      id: "Oligodendrocytes",
      negEntropy: 0,
    },
    {
      id: PERIPHERAL_SENSORY,
      negEntropy: 0,
    },
    {
      id: SPINAL_CORD,
      negEntropy: 0,
    },
    {
      id: "Sympathetic",
      negEntropy: 0,
    },
    {
      id: TELENCEPHALON_INTERNEURONS,
      negEntropy: 0,
    },
    {
      id: TELENCEPHALON_PROJECTING,
      negEntropy: 0,
    },
    {
      id: "Vascular",
      negEntropy: 0,
    }, // Di- and mesencephalon children
    {
      id: "Excitatory",
      negEntropy: 0.7,
      parentId: DI_AND_MESENCEPHALON,
    },
    {
      id: "Inhibitory",
      negEntropy: 0.7,
      parentId: DI_AND_MESENCEPHALON,
    }, // Excitatory cells
    {
      id: "CR",
      parentId: "Excitatory",
      proportionalExpression: 0.9,
      relativeExpression: 0.85,
    },
    {
      id: "DECHO2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU4",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU5",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU4",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU5",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU6",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU7",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU8",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU9",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU10",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU11",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
  ],
  name: "Nhlh2",
};

export const GM5741 = {
  data: [
    {
      id: "Astroependymal",
      negEntropy: 0,
    },
    {
      id: "Cerebellum",
      negEntropy: 0,
    },
    {
      id: CHOLINERGIC_MONOAMINERGIC,
      negEntropy: 0,
    },
    {
      id: DI_AND_MESENCEPHALON,
      negEntropy: 0.8,
    },
    {
      id: "Enteric",
      negEntropy: 0,
    },
    {
      id: "Hindbrain",
      negEntropy: 0,
    },
    {
      id: IMMATURE_NEURAL,
      negEntropy: 0,
    },
    {
      id: "Immune",
      negEntropy: 0,
    },
    {
      id: NEURAL_CREST_LIKE_GLIA,
      negEntropy: 0,
    },
    {
      id: "Oligodendrocytes",
      negEntropy: 0,
    },
    {
      id: PERIPHERAL_SENSORY,
      negEntropy: 0,
    },
    {
      id: SPINAL_CORD,
      negEntropy: 0,
    },
    {
      id: "Sympathetic",
      negEntropy: 0,
    },
    {
      id: TELENCEPHALON_INTERNEURONS,
      negEntropy: 0,
    },
    {
      id: TELENCEPHALON_PROJECTING,
      negEntropy: 0,
    },
    {
      id: "Vascular",
      negEntropy: 0,
    }, // Di- and mesencephalon children
    {
      id: "Excitatory",
      negEntropy: 0.8,
      parentId: DI_AND_MESENCEPHALON,
    },
    {
      id: "Inhibitory",
      negEntropy: 0.4,
      parentId: DI_AND_MESENCEPHALON,
    },
    // Excitatory cells
    {
      id: "CR",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DECHO2",
      parentId: "Excitatory",
      proportionalExpression: 0.65,
      relativeExpression: 0.8,
    },
    {
      id: "DEGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU4",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU5",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU4",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU5",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU6",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU7",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU8",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU9",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU10",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU11",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
  ],
  name: "Gm5741",
};

export const GNG8 = {
  data: [
    {
      id: "Astroependymal",
      negEntropy: 0,
    },
    {
      id: "Cerebellum",
      negEntropy: 0,
    },
    {
      id: CHOLINERGIC_MONOAMINERGIC,
      negEntropy: 0,
    },
    {
      id: DI_AND_MESENCEPHALON,
      negEntropy: 0.8,
    },
    {
      id: "Enteric",
      negEntropy: 0,
    },
    {
      id: "Hindbrain",
      negEntropy: 0,
    },
    {
      id: IMMATURE_NEURAL,
      negEntropy: 0,
    },
    {
      id: "Immune",
      negEntropy: 0,
    },
    {
      id: NEURAL_CREST_LIKE_GLIA,
      negEntropy: 0,
    },
    {
      id: "Oligodendrocytes",
      negEntropy: 0,
    },
    {
      id: PERIPHERAL_SENSORY,
      negEntropy: 0,
    },
    {
      id: SPINAL_CORD,
      negEntropy: 0,
    },
    {
      id: "Sympathetic",
      negEntropy: 0,
    },
    {
      id: TELENCEPHALON_INTERNEURONS,
      negEntropy: 0,
    },
    {
      id: TELENCEPHALON_PROJECTING,
      negEntropy: 0,
    },
    {
      id: "Vascular",
      negEntropy: 0,
    }, // Di- and mesencephalon children
    {
      id: "Excitatory",
      negEntropy: 0.7,
      parentId: DI_AND_MESENCEPHALON,
    },
    {
      id: "Inhibitory",
      negEntropy: 0.1,
      parentId: DI_AND_MESENCEPHALON,
    }, // Excitatory cells
    {
      id: "CR",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DECHO2",
      parentId: "Excitatory",
      proportionalExpression: 0.9,
      relativeExpression: 0.7,
    },
    {
      id: "DEGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU4",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU5",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU4",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU5",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU6",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU7",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU8",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU9",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU10",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU11",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
  ],
  name: "Gng8",
};

export const LRRC55 = {
  data: [
    {
      id: "Astroependymal",
      negEntropy: 0,
    },
    {
      id: "Cerebellum",
      negEntropy: 0,
    },
    {
      id: CHOLINERGIC_MONOAMINERGIC,
      negEntropy: 0,
    },
    {
      id: DI_AND_MESENCEPHALON,
      negEntropy: 0.8,
    },
    {
      id: "Enteric",
      negEntropy: 0,
    },
    {
      id: "Hindbrain",
      negEntropy: 0,
    },
    {
      id: IMMATURE_NEURAL,
      negEntropy: 0,
    },
    {
      id: "Immune",
      negEntropy: 0,
    },
    {
      id: NEURAL_CREST_LIKE_GLIA,
      negEntropy: 0,
    },
    {
      id: "Oligodendrocytes",
      negEntropy: 0,
    },
    {
      id: PERIPHERAL_SENSORY,
      negEntropy: 0,
    },
    {
      id: SPINAL_CORD,
      negEntropy: 0,
    },
    {
      id: "Sympathetic",
      negEntropy: 0,
    },
    {
      id: TELENCEPHALON_INTERNEURONS,
      negEntropy: 0,
    },
    {
      id: TELENCEPHALON_PROJECTING,
      negEntropy: 0,
    },
    {
      id: "Vascular",
      negEntropy: 0,
    }, // Di- and mesencephalon children
    {
      id: "Excitatory",
      negEntropy: 0.7,
      parentId: DI_AND_MESENCEPHALON,
    },
    {
      id: "Inhibitory",
      negEntropy: 0.2,
      parentId: DI_AND_MESENCEPHALON,
    },
    // Excitatory cells
    {
      id: "CR",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DECHO2",
      parentId: "Excitatory",
      proportionalExpression: 0.9,
      relativeExpression: 0.7,
    },
    {
      id: "DEGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU4",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU5",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU1",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU2",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU3",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU4",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU5",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU6",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU7",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU8",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU9",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU10",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU11",
      parentId: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
  ],
  name: "Lrrc55",
};

export const GENES = [GM5741, GNG8, LRRC55, NDNF, NHLH2];
