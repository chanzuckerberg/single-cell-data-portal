import { Dialog } from "@czi-sds/components";
import { FC } from "react";

interface Props {
  open: boolean;
  onClose: () => void;
}

const BYODModal: FC<Props> = ({ open, onClose }) => {
  return (
    <Dialog open={open} onClose={onClose} title="Explore Your Data">
      <div>Copy pending approval</div>
    </Dialog>
  );
};

export default BYODModal;
