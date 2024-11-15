import { SvgIcon } from "@mui/material";
import { Props } from "src/views/Collection/components/CollectionActions/components/MoreDropdown/components/Menu/components/IconSort/types";

export default function IconSort({ className }: Props): JSX.Element {
  return (
    <SvgIcon className={className} color="inherit">
      <svg
        fill="none"
        height="12"
        viewBox="0 0 12 12"
        width="12"
        xmlns="http://www.w3.org/2000/svg"
      >
        <path
          d="M3.70711 0.292893C3.31658 -0.0976311 2.68342 -0.0976311 2.29289 0.292893L0.292893 2.29289C-0.0976311 2.68342 -0.0976311 3.31658 0.292893 3.70711C0.683417 4.09763 1.31658 4.09763 1.70711 3.70711L2 3.41421V10C2 10.5523 2.44772 11 3 11C3.55228 11 4 10.5523 4 10V3.41421L4.29289 3.70711C4.68342 4.09763 5.31658 4.09763 5.70711 3.70711C6.09763 3.31658 6.09763 2.68342 5.70711 2.29289L3.70711 0.292893Z"
          fill="currentColor"
        />
        <path
          d="M11.7071 9.70711L9.70711 11.7071C9.31658 12.0976 8.68342 12.0976 8.29289 11.7071L6.29289 9.70711C5.90237 9.31658 5.90237 8.68342 6.29289 8.29289C6.68342 7.90237 7.31658 7.90237 7.70711 8.29289L8 8.58579V2C8 1.44771 8.44772 1 9 1C9.55228 1 10 1.44771 10 2V8.58579L10.2929 8.29289C10.6834 7.90237 11.3166 7.90237 11.7071 8.29289C12.0976 8.68342 12.0976 9.31658 11.7071 9.70711Z"
          fill="currentColor"
        />
      </svg>
    </SvgIcon>
  );
}
