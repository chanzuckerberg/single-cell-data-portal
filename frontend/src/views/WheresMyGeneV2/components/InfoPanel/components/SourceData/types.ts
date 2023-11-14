export interface Collection {
  name: string;
  url: string;
  datasets: { id: string; label: string }[];
  total_count: number;
}

export interface Collections {
  [name: string]: Collection;
}
