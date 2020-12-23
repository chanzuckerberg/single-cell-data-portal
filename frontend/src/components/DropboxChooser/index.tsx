import loadable from "@loadable/component";
import React, { FC, useState } from "react";
import { Dataset } from "src/common/entities";

declare global {
  interface Window {
    Dropbox: {
      isBrowserSupported: () => boolean;
      choose: (options: unknown) => void;
    };
  }
}

const AsyncContent = loadable(
  () =>
    /*webpackChunkName: 'DropboxChooser/Content' */ import(
      "src/components/DropboxChooser/components/Content"
    )
);

export interface DropboxFile {
  link: string;
  name: string;
}

export type UploadingFile = DropboxFile & Partial<Dataset>;

const DROPBOX_OPTIONS = {
  extensions: [".h5ad"],
  sizeLimit: 30 * 10 ** 30, // 30GB
};

export interface Props {
  onUploadFile: (newFile: UploadingFile) => void;
}

const DropboxChooser: FC<Props> = ({ children, onUploadFile }) => {
  const [isContentShown, setIsContentShown] = useState(false);

  const handleMouseOver = () => {
    setIsContentShown(true);
  };

  const handleClick = () => {
    if (!window.Dropbox.isBrowserSupported()) {
      alert("Sorry, your browser does not support Dropbox upload :(");
    }

    const options = {
      ...DROPBOX_OPTIONS,
      success(files: DropboxFile[] = []) {
        onUploadFile({ id: "", link: files[0].link, name: files[0].name });
      },
    };

    window.Dropbox.choose(options);
  };

  return (
    <>
      <div onMouseOver={handleMouseOver} onClick={handleClick}>
        {children}
      </div>
      {isContentShown && <AsyncContent />}
    </>
  );
};

export default DropboxChooser;
