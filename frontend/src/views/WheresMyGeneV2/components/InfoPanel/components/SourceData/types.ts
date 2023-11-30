export interface Collection {
  name: string;
  url: string;
  datasets: { id: string; label: string }[];
}

export interface Collections {
  [name: string]: Collection;
}
