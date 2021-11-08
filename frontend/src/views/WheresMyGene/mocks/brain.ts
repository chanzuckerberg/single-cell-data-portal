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
  },
  {
    id: "Cerebellum",
    name: "Cerebellum",
  },
  {
    id: CHOLINERGIC_MONOAMINERGIC,
    name: CHOLINERGIC_MONOAMINERGIC,
  },
  {
    id: DI_AND_MESENCEPHALON,
    name: DI_AND_MESENCEPHALON,
  },
  {
    id: "Enteric",
    name: "Enteric",
  },
  {
    id: "Hindbrain",
    name: "Hindbrain",
  },
  {
    id: IMMATURE_NEURAL,
    name: IMMATURE_NEURAL,
  },
  {
    id: "Immune",
    name: "Immune",
  },
  {
    id: NEURAL_CREST_LIKE_GLIA,
    name: NEURAL_CREST_LIKE_GLIA,
  },
  {
    id: "Oligodendrocytes",
    name: "Oligodendrocytes",
  },
  {
    id: PERIPHERAL_SENSORY,
    name: PERIPHERAL_SENSORY,
  },
  {
    id: SPINAL_CORD,
    name: SPINAL_CORD,
  },
  {
    id: "Sympathetic",
    name: "Sympathetic",
  },
  {
    id: TELENCEPHALON_INTERNEURONS,
    name: TELENCEPHALON_INTERNEURONS,
  },
  {
    id: TELENCEPHALON_PROJECTING,
    name: TELENCEPHALON_PROJECTING,
  },
  {
    id: "Vascular",
    name: "Vascular",
  },
  // Di- and mesencephalon children
  {
    id: "Excitatory",
    name: "Excitatory",
  },
  {
    id: "Inhibitory",
    name: "Inhibitory",
  },
  // Excitatory cells
  {
    id: "CR",
    name: "CR",
  },
  {
    id: "DECHO2",
    name: "DECHO2",
  },
  {
    id: "DEGLU1",
    name: "DEGLU1",
  },
  {
    id: "DEGLU2",
    name: "DEGLU2",
  },
  {
    id: "DEGLU3",
    name: "DEGLU3",
  },
  {
    id: "DEGLU4",
    name: "DEGLU4",
  },
  {
    id: "DEGLU5",
    name: "DEGLU5",
  },
  {
    id: "HBGLU1",
    name: "HBGLU1",
  },
  {
    id: "HBGLU2",
    name: "HBGLU2",
  },
  {
    id: "HBGLU3",
    name: "HBGLU3",
  },
  {
    id: "MEGLU1",
    name: "MEGLU1",
  },
  {
    id: "MEGLU2",
    name: "MEGLU2",
  },
  {
    id: "MEGLU3",
    name: "MEGLU3",
  },
  {
    id: "MEGLU4",
    name: "MEGLU4",
  },
  {
    id: "MEGLU5",
    name: "MEGLU5",
  },
  {
    id: "MEGLU6",
    name: "MEGLU6",
  },
  {
    id: "MEGLU7",
    name: "MEGLU7",
  },
  {
    id: "MEGLU8",
    name: "MEGLU8",
  },
  {
    id: "MEGLU9",
    name: "MEGLU9",
  },
  {
    id: "MEGLU10",
    name: "MEGLU10",
  },
  {
    id: "MEGLU11",
    name: "MEGLU11",
  },
];

// GENES
export const NDNF = {
  data: [
    {
      id: "Astroependymal",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Cerebellum",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: CHOLINERGIC_MONOAMINERGIC,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: DI_AND_MESENCEPHALON,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Enteric",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Hindbrain",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: IMMATURE_NEURAL,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Immune",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: NEURAL_CREST_LIKE_GLIA,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Oligodendrocytes",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: PERIPHERAL_SENSORY,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: SPINAL_CORD,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Sympathetic",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: TELENCEPHALON_INTERNEURONS,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: TELENCEPHALON_PROJECTING,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Vascular",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    // Di- and mesencephalon children
    {
      id: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Inhibitory",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    // Excitatory cells
    {
      id: "CR",
      proportionalExpression: 0.9,
      relativeExpression: 0.85,
    },
    {
      id: "DECHO2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU4",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU5",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU4",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU5",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU6",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU7",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU8",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU9",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU10",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU11",
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
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Cerebellum",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: CHOLINERGIC_MONOAMINERGIC,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: DI_AND_MESENCEPHALON,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Enteric",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Hindbrain",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: IMMATURE_NEURAL,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Immune",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: NEURAL_CREST_LIKE_GLIA,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Oligodendrocytes",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: PERIPHERAL_SENSORY,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: SPINAL_CORD,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Sympathetic",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: TELENCEPHALON_INTERNEURONS,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: TELENCEPHALON_PROJECTING,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Vascular",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    }, // Di- and mesencephalon children
    {
      id: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Inhibitory",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    }, // Excitatory cells
    {
      id: "CR",
      proportionalExpression: 0.9,
      relativeExpression: 0.85,
    },
    {
      id: "DECHO2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU4",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU5",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU4",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU5",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU6",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU7",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU8",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU9",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU10",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU11",
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
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Cerebellum",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: CHOLINERGIC_MONOAMINERGIC,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: DI_AND_MESENCEPHALON,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Enteric",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Hindbrain",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: IMMATURE_NEURAL,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Immune",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: NEURAL_CREST_LIKE_GLIA,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Oligodendrocytes",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: PERIPHERAL_SENSORY,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: SPINAL_CORD,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Sympathetic",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: TELENCEPHALON_INTERNEURONS,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: TELENCEPHALON_PROJECTING,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Vascular",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    }, // Di- and mesencephalon children
    {
      id: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Inhibitory",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    // Excitatory cells
    {
      id: "CR",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DECHO2",
      proportionalExpression: 0.65,
      relativeExpression: 0.8,
    },
    {
      id: "DEGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU4",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU5",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU4",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU5",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU6",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU7",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU8",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU9",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU10",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU11",
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
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Cerebellum",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: CHOLINERGIC_MONOAMINERGIC,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: DI_AND_MESENCEPHALON,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Enteric",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Hindbrain",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: IMMATURE_NEURAL,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Immune",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: NEURAL_CREST_LIKE_GLIA,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Oligodendrocytes",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: PERIPHERAL_SENSORY,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: SPINAL_CORD,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Sympathetic",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: TELENCEPHALON_INTERNEURONS,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: TELENCEPHALON_PROJECTING,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Vascular",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    }, // Di- and mesencephalon children
    {
      id: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Inhibitory",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    }, // Excitatory cells
    {
      id: "CR",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DECHO2",
      proportionalExpression: 0.9,
      relativeExpression: 0.7,
    },
    {
      id: "DEGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU4",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU5",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU4",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU5",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU6",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU7",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU8",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU9",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU10",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU11",
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
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Cerebellum",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: CHOLINERGIC_MONOAMINERGIC,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: DI_AND_MESENCEPHALON,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Enteric",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Hindbrain",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: IMMATURE_NEURAL,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Immune",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: NEURAL_CREST_LIKE_GLIA,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Oligodendrocytes",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: PERIPHERAL_SENSORY,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: SPINAL_CORD,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Sympathetic",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: TELENCEPHALON_INTERNEURONS,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: TELENCEPHALON_PROJECTING,
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Vascular",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    }, // Di- and mesencephalon children
    {
      id: "Excitatory",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    {
      id: "Inhibitory",
      proportionalExpression: 0.35,
      relativeExpression: 0.12,
    },
    // Excitatory cells
    {
      id: "CR",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DECHO2",
      proportionalExpression: 0.9,
      relativeExpression: 0.7,
    },
    {
      id: "DEGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU4",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "DEGLU5",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "HBGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU1",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU2",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU3",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU4",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU5",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU6",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU7",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU8",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU9",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU10",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
    {
      id: "MEGLU11",
      proportionalExpression: 0.35,
      relativeExpression: 0.2,
    },
  ],
  name: "Lrrc55",
};

export const GENES = [GM5741, GNG8, LRRC55, NDNF, NHLH2];
